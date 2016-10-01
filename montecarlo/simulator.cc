/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "simulator.h"

namespace mc {

Simulator::Simulator(input::Parameters& parms) : LatticeGraph(parms)
{
  Model::construct(LatticeGraph::lattice(), parms);

  // system state
  state.resize(num_sites());
  for (auto s=sites_begin(); s!=sites_end(); ++s) {
    unsigned stype = site_type(s);
    state_idx max_idx = sitebasis_dimension(stype)-1;
    if (max_idx < 1) {
      throw std::range_error("Simulator::Simulator: site basis dimension must be >= 2");
    }
    state[site(s)] = SiteState(stype, 0, max_idx);
  }

  // Random number generator
  int seed_type = parms.set_value("rng_seed", 0);
  rng.seed(seed_type);
  // lattice site generator 
  rng.set_site_generator(0, num_sites()-1);
  // site types & a RNG of state indices for each site type
  for (const unsigned& site_type : site_types()) {
    unsigned max_idx = sitebasis_dimension(site_type)-1;
    rng.add_state_generator(site_type, 0, max_idx);
  }

  // MC parameters
  measure_samples = parms.set_value("measure_samples", 0);
  warmup = parms.set_value("warmup", 0);
  min_interval = parms.set_value("min_interval", 0);

  // observables
  observables.init(parms, *this, print_copyright);
  //init_observables(parms);
}

void Simulator::start(input::Parameters& parms) 
{
  // update model parameters
  Model::update_parameters(parms);

  // set variable objects
  T = parms.set_value("T", 1.0);
  kB = parms.set_value("kB", 1.0);
  beta = 1.0/(kB*T);

  // values of exp(-beta*E) for a bond for Metropolis algorithm
  init_boltzmann_table();

  // observables
  observables.reset();

  /*-----------------Simulation START-----------------*/
  // initial state
  init_state_random();
  // warmup runs
  for (int i=0; i<warmup; ++i) update_state_metropolis();

  // measurement
  int count = 0;
  int num_measurement = 0;
  while (num_measurement != measure_samples) {
    update_state_metropolis();
    if (count++ == min_interval) {
      // make measurements
      count = 0;
      num_measurement++;
      do_measurements();
    }
  }
  // output
  observables.print(T);
  /*-----------------Simulation END-----------------*/
}

// measurements
inline void Simulator::do_measurements(void)
{
  // energy
  if (observables.energy_terms()) {
    observables.energy_terms() << get_energy();
  }

  // magnetization
  if (observables.magn() || observables.magn_sq()) {
    double x = magnetization();
    if (observables.magn()) observables.magn() << x;
    if (observables.magn_sq()) observables.magn_sq() << x*x;
  }

  // Potts magnetization
  if (observables.potts_magn() || observables.potts_magn_sq()) {
    double x = potts_magnetization();
    if (observables.potts_magn()) observables.potts_magn() << x;
    if (observables.potts_magn_sq()) observables.potts_magn_sq() << x*x;
  }

  // strain
  if (observables.strain() || observables.strain_sq()) {
    double x = strain();
    if (observables.strain()) observables.strain() << x;
    if (observables.strain_sq()) observables.strain_sq() << x * x;
  }
}

void Simulator::init_state_random(void)
{
  // random initial state
  for (auto& s : state) {
    state_idx idx = rng.random_idx(s.type());
    s.reset_idx(idx);
  }
} 

void Simulator::update_state_metropolis(void)
{
  out_bond_iterator ob, ob_end;
  in_bond_iterator ib, ib_end;
  for (unsigned i=0; i<num_sites(); ++i) {
    // suggest new state and random site
    unsigned site = rng.random_site();
    unsigned curr_idx = state[site].idx();
    unsigned new_idx = curr_idx;
    while (new_idx==curr_idx) new_idx = rng.random_idx(state[site].type());

    double w_old = 1.0;
    double w_new = 1.0;
    for (std::tie(ob,ob_end)=out_bonds(site); ob!=ob_end; ++ob) {
      w_old *= boltzmann_table[bond_type(ob)](curr_idx, state[target(ob)].idx());
      w_new *= boltzmann_table[bond_type(ob)](new_idx, state[target(ob)].idx());
    }
    for (std::tie(ib,ib_end)=in_bonds(site); ib!=ib_end; ++ib) {
      w_old *= boltzmann_table[bond_type(ib)](state[source(ib)].idx(), curr_idx);
      w_new *= boltzmann_table[bond_type(ib)](state[source(ib)].idx(), new_idx);
    }
    double proby = w_new/w_old;

    // acceptance
    if (proby > rng.random_real()) {
      state[site].reset_idx(new_idx);
    }
  }
}


// exp(-beta*E) factors for different bond types
void Simulator::init_boltzmann_table(void)
{
  //boltzmann_table.clear();
  for (auto& mat : boltzmann_table) mat = Eigen::MatrixXd();
  // get bond types
  for (unsigned i=0; i<lattice().num_unitcell_bonds(); ++i) {
    lattice::Bond b = lattice().unitcell_bond(i);
    lattice::Site src = lattice().unitcell_site(b.src_id());
    lattice::Site tgt = lattice().unitcell_site(b.tgt_id());
    unsigned src_dim = sitebasis_dimension(src.type());
    unsigned tgt_dim = sitebasis_dimension(tgt.type());
    Eigen::MatrixXd exp_betaE(src_dim, tgt_dim);
    // energy for this two sites
    for (unsigned i=0; i<src_dim; ++i) {
      // siteterm for src site
      double E_src = 0.0;
      for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
        double m = sterm->matrix_element(src.type(), i);
        E_src += m * sterm->coupling(src.type());
      }
      for (unsigned j=0; j<tgt_dim; ++j) {
        // siteterm for src site
        double E_tgt = 0.0;
        for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
          double m = sterm->matrix_element(tgt.type(), j);
          E_tgt += m * sterm->coupling(tgt.type());
        }
        // bond energy
        double E_bond = 0.0;
        for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
          double m = bterm->matrix_element(b.type(), i, j);
          E_bond += m * bterm->coupling(b.type());
        }
        exp_betaE(i,j) = std::exp(-beta*(E_src+E_tgt+E_bond));
      }
    }
    // matrix built for this bond type
    //boltzmann_table.insert({b.type(), exp_betaE});
    boltzmann_table[b.type()] = exp_betaE;
  }
}

void Simulator::init_observables(const input::Parameters& parms)
{
  /*need_energy = parms.set_value("energy", false);
  need_magn = parms.set_value("magnetization", false);
  if (need_magn) {
    magn_op.init(basis(), "S(i)");
    magn_data.init("Magnetization");
    magn_out.open("res_magnetization.txt");
    magn_out << std::left;
    magn_out << std::setw(12) << "#T";
    magn_out << std::setw(12) << "Magn";
    magn_out << std::setw(12) << "err";
    magn_out << std::setw(12) << "samples";
    magn_out << std::setw(12) << "conv";
    magn_out << std::setw(12) << "tau";
    magn_out << std::endl;
  }
  if (need_energy) {
    energy_size = num_siteterms() + num_bondterms();
    energy_data.init("Energy", energy_size);
    energy_out.open("res_energy.txt");
  }
  */
}


/*void Simulator::init(void)
{
  // random initial state
  init_state_random();

  // energy
  //si_end = std::vector<SiteTerm>::end();
  for (auto s=sites_begin(); s!=sites_end(); ++s) {
    state_idx idx = state[site(s)].idx();
    unsigned type = state[site(s)].type();
    for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
      double m = sterm->matrix_element(type, idx);
      //std::cout << *s << " " << m << std::endl;
    }
  }

  // Bond Terms
  for (auto b=bonds_begin(); b!=bonds_end(); ++b) {
    unsigned type = bond_type(b);
    site_descriptor src = source(b);
    site_descriptor tgt = target(b);
    state_idx src_idx = state[site(src)].idx();
    state_idx tgt_idx = state[site(tgt)].idx();
    for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
      double m = bterm->matrix_element(type, src_idx, tgt_idx);
      //std::cout << "bond: " << " " << m << std::endl;
    }
  }

  //for (unsigned i=0; i<state.size(); ++i) {
    //unsigned site_type = vertex_type(i);
    //state_idx idx = state[i].idx();
   // for (std::tie(st,st_end)=site_terms(); st!=st_end; ++st) {
      //double m = st->matrix_element(site_type, idx);
      //std::cout << i << " " << m << std::endl;
   // }
    //}
    //std::cout << vertex_id(vi) << std::endl;
    //matrixelem[site_type] += siteterm[i].apply_operator(site_type, state_idx);
    // evaluate(siteterm site_type, stat_idx)
    //}
  //}
}*/ 
void Simulator::print_copyright(std::ostream& os) 
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: Monte Carlo Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n"; 
  os << "#" << std::string(72,'-') << "\n";
}


} // end namespace basis
