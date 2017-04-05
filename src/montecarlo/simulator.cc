/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-03 18:28:22
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "simulator.h"

namespace mc {

Simulator::Simulator(input::Parameters& parms) 
  : LatticeGraph(parms), Model(parms, LatticeGraph::lattice())
{
  //Model::construct(LatticeGraph::lattice(), parms);

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

  // Boltzmann constant
  kB = parms.set_value("kB", 1.0);

  // observables
  observables.init(parms, *this, print_copyright);
  // extra consideration for energy & magnetization
  need_energy = false;
  if (observables.energy_terms() || observables.energy() || observables.energy_sq()) {
    energy_terms.resize(num_total_terms());
    need_energy = true;
  }
  if (observables.energy_terms_sq()) {
    energy_terms_sq.resize(num_total_terms());
    need_energy = true;
  }
  need_magn = false;
  if (observables.magn() || observables.magn_sq() || 
    observables.potts_magn() || observables.potts_magn_sq()) {
    need_magn = true;
  }
}

int Simulator::start(input::Parameters& parms) 
{
  // this function must be override
  throw std::logic_error("Simulator::start: You must override this function");

  //Model::update_parameters(parms);
  //// set variable parameters
  //T = parms.set_value("T", 1.0);
  //kB = parms.set_value("kB", 1.0);
  //beta = 1.0/(kB*T);
  //// values of exp(-beta*E) for bond energy (E) for Metropolis algorithm
  //init_boltzmann_table();
  //// observables
  //observables.reset();
  ///*-----------------Simulation START-----------------*/
  //// initial state
  //init_state_random();
  //// warmup runs
  //for (int i=0; i<warmup; ++i) update_state_metropolis();
  //// measurement
  //int count = 0;
  //int num_measurement = 0;
  //while (num_measurement != measure_samples) {
  //  update_state_metropolis();
  //  if (count++ == min_interval) {
  //    // make measurements
  //    count = 0;
  //    num_measurement++;
  //    do_measurements();
  //  }
  //}
  //// output
  //observables.print(T);
  /*-----------------Simulation END-----------------*/
  return 0;
}

void Simulator::update_parameters(input::Parameters& parms) {
  // update model parameters
  Model::update_parameters(parms);

  // set variable parameters
  T = parms.set_value("T", 1.0);
  beta = 1.0/(kB*T);

  // values of exp(-beta*E) for bond energy (E) for Metropolis algorithm
  init_boltzmann_table();

  std::cout << "Simulator:: hack in calculating Potts magnetization\n";
}

/*----------------------Energy-----------------------*/
inline mc::VectorData Simulator::get_energy(void) const
{
  unsigned num_terms = num_total_terms();
  mc::VectorData energy(num_terms);
  for (unsigned i=0; i<num_terms; ++i) energy(i) = 0.0;
  // bond energies
  for (auto b=bonds_begin(); b!=bonds_end(); ++b) {
    unsigned type = bond_type(b);
    auto src = source(b);
    auto tgt = target(b);
    state_idx src_idx = state[site(src)].idx();
    state_idx tgt_idx = state[site(tgt)].idx();
    unsigned term = 0;
    for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
      //std::cout << src << " - " << tgt << "\n";
      //std::cout << strain_op.apply(state[src]) << " - ";
      //std::cout << strain_op.apply(state[tgt]) << "\n";
      double m = bterm->matrix_element(type, src_idx, tgt_idx);
      double c = bterm->coupling(type);
      //std::cout << "bond: " << " " << m << std::endl;
      energy[term++] += m * c;
      //std::cout << "energy " << m*c << "\n";
    }
    //getchar();
  }
  // site energies
  for (const auto& s : state) {
    unsigned term = num_bondterms();
    for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
      double m = sterm->matrix_element(s.type(), s.idx());
      double c = sterm->coupling(s.type());
      energy[term++] += m * c;
    }
  }
  // energy per site
  return energy/num_sites();
}

/*----------------------Ising magnetization-----------------------*/
inline double Simulator::get_magnetization(void)
{
  double ms = 0;
  for (const auto& s : state) {
    ms += magn_op.apply(s);
  }
  return std::abs(ms)/num_sites();
}

/*----------------------Potts magnetization-----------------------*/
inline double Simulator::get_potts_magnetization(void)
{
  unsigned num_sitetypes = basis().size();
  std::vector<std::multiset<int> > lattice_spins(num_sitetypes);
  std::vector<std::set<int> > site_spins(num_sitetypes);

  // spin values for a sites
  // HACK
  //for (unsigned i=0; i<num_sitetypes; ++i) {
  for (unsigned i=1; i<num_sitetypes; ++i) {
    for (unsigned j=0; j<sitebasis_dimension(i); ++j) {
      SiteBasisState site_state(i, j, sitebasis_dimension(i)-1); 
      int spin = std::nearbyint(magn_op.apply(site_state));
      //std::cout << "i, spin = " << i << ", " << spin << "\n";
      site_spins[i].insert(spin);
    }
  }
  //for (const auto& s : site_spins[0]) std::cout << s << "\n";

  // set of spin values for all the sites in the lattice
  for (const auto& s : state) {
    int spin = std::nearbyint(magn_op.apply(s));
    lattice_spins[s.type()].insert(spin);
  }

  // magnetization
  double ms = 0;
  // HACK
  //for (unsigned i=0; i<num_sitetypes; ++i) {
  for (unsigned i=1; i<num_sitetypes; ++i) {
    int q = site_spins[i].size();
    int Nmax = 0;
    for (const auto& s : site_spins[i]) {
      int count = lattice_spins[i].count(s);
      //std::cout << "s, count = " << s << ", " << count << "\n";
      if (Nmax < count) Nmax = count;
    }
    //std::cout << "Nmax, Nmag = " << Nmax << ", " << lattice_spins[i].size() << "\n";
    ms += static_cast<double>(q*Nmax - lattice_spins[i].size())/(q-1);
  }
  //std::cout << "M = " << ms/num_sites() << "\n";
  //abort();
  return ms/num_sites();
}

/*----------------------Strain-----------------------*/
inline double Simulator::get_strain(void)
{
  double ms = 0;
  for (const auto& s : state) {
    ms += strain_op.apply(s);
    //std::cout << "strain = " << strain_op.apply(s); getchar();
  }
  return std::abs(ms)/num_sites();
  //return ms/num_sites();
}


/*------------------random initial state----------------*/
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
    unsigned site_type = state[site].type();
    unsigned curr_idx = state[site].idx();
    unsigned new_idx = curr_idx;
    while (new_idx==curr_idx) new_idx = rng.random_idx(site_type);

    double w_old = 1.0;
    double w_new = 1.0;
    //std::cout << "site, sigma = " << site << " " << strain_op.apply(state[site]) << "\n";
    for (std::tie(ob,ob_end)=out_bonds(site); ob!=ob_end; ++ob) {
      //std::cout << boltzmann_table[bond_type(ob)].cols() << "\n";
      //std::cout << curr_idx << state[target(ob)].idx() << "\n";
      //std::cout << "sigma = " << strain_op.apply(state[target(ob)]) << "\n";
      w_old *= boltzmann_table[bond_type(ob)](curr_idx, state[target(ob)].idx());
      w_new *= boltzmann_table[bond_type(ob)](new_idx, state[target(ob)].idx());
    }
    for (std::tie(ib,ib_end)=in_bonds(site); ib!=ib_end; ++ib) {
      //std::cout << bond_type(ib) << " " << state[source(ib)].idx() << " " << curr_idx << "\n";
      //std::cout << boltzmann_table[bond_type(ib)].rows() << "\n";
      //std::cout << "sigma = " << strain_op.apply(state[source(ib)]) << "\n";
      w_old *= boltzmann_table[bond_type(ib)](state[source(ib)].idx(), curr_idx);
      w_new *= boltzmann_table[bond_type(ib)](state[source(ib)].idx(), new_idx);
    }
    double proby = w_new/w_old;
    //std::cout << "site = " << site << "\n";
    //std::cout << "proby = " << proby << "\n";
    //getchar();

    // acceptance
    if (proby > rng.random_real()) {
      state[site].reset_idx(new_idx);
    }
  }
}

// exp(-beta*E) factors for different bond types
void Simulator::init_boltzmann_table(void)
{
  for (auto& mat : boltzmann_table) mat = Eigen::MatrixXd();
  // get bond types
  //for (unsigned i=0; i<lattice().num_unitcell_bonds(); ++i) {
  //  lattice::Bond b = lattice().unitcell_bond(i);
  //  lattice::Site src = lattice().unitcell_site(b.src_id());
  //  lattice::Site tgt = lattice().unitcell_site(b.tgt_id());
  unsigned btype, src_type, tgt_type;
  for (const auto& elem : Model::bond_sites_map()) {
    btype = elem.first;
    std::tie(src_type, tgt_type) = elem.second;
    unsigned src_dim = sitebasis_dimension(src_type);
    unsigned tgt_dim = sitebasis_dimension(tgt_type);
    Eigen::MatrixXd exp_betaE(src_dim, tgt_dim);
    // energy for this two sites
    for (unsigned i=0; i<src_dim; ++i) {
      // siteterm for src site
      double E_src = 0.0;
      for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
        double m = sterm->matrix_element(src_type, i);
        E_src += m * sterm->coupling(src_type);
      }
      for (unsigned j=0; j<tgt_dim; ++j) {
        // siteterm for src site
        double E_tgt = 0.0;
        for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
          double m = sterm->matrix_element(tgt_type, j);
          E_tgt += m * sterm->coupling(tgt_type);
        }
        // bond energy
        double E_bond = 0.0;
        for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
          double m = bterm->matrix_element(btype, i, j);
          E_bond += m * bterm->coupling(btype);
        }
        exp_betaE(i,j) = std::exp(-beta*(E_src+E_tgt+E_bond));
      }
    }
    // matrix built for this bond type
    boltzmann_table[btype] = exp_betaE;
  }
}

void Simulator::print_copyright(std::ostream& os) 
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: Monte Carlo Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n"; 
  os << "#" << std::string(72,'-') << "\n";
}

// measurements
/*inline void Simulator::do_measurements(void)
{
  // energy
  if (need_energy) {
    get_energy(energy_terms);
    if (observables.energy_terms()) {
      observables.energy_terms() << energy_terms;
    }
    if (observables.energy()) {
      double e = energy_terms.sum();
      observables.energy() << e;
    }
    if (observables.energy_sq()) {
      double e = energy_terms.sum();
      observables.energy_sq() << e*e;
    }
    if (observables.energy_terms_sq()) {
      observables.energy_terms_sq() << energy_terms.square();
    }
  }

  // magnetization
  if (observables.magn() || observables.magn_sq()) {
    double x = get_magnetization();
    if (observables.magn()) observables.magn() << x;
    if (observables.magn_sq()) observables.magn_sq() << x*x;
  }

  // Potts magnetization
  if (observables.potts_magn() || observables.potts_magn_sq()) {
    double x = get_potts_magnetization();
    if (observables.potts_magn()) observables.potts_magn() << x;
    if (observables.potts_magn_sq()) observables.potts_magn_sq() << x*x;
  }

  // strain
  if (observables.strain() || observables.strain_sq()) {
    double x = get_strain();
    if (observables.strain()) observables.strain() << x;
    if (observables.strain_sq()) observables.strain_sq() << x * x;
  }
}*/


} // end namespace basis
