/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "mc_ni2mnx.h"

MCSimulation::MCSimulation(input::Parameters& parms) : Simulator(parms) 
{
  // observables as functions of what?
  observables.as_function_of("T^*");

  // define your observable operators (other than energy)
  if (need_magn) {
    Simulator::set_magn_op("S(i)");
  }
  if (observables.strain() || observables.strain_sq()) {
    Simulator::set_strain_op("sigma(i)");
  }

  // Mn2 doping
  double mn2_percent = parms.set_value("Mn2_doping_percent", 0.0);
  if (mn2_percent<0.0 || mn2_percent>100.0) 
    throw std::out_of_range("Invalid value of 'Mn2_doping_percent'");
  // choose the 40% of sites randomly to dope with Mn2 
  unsigned num_doped_sites = std::nearbyint(0.01 * mn2_percent * num_sites());
  if (num_doped_sites > 0) {
    std::vector<unsigned> all_sites(num_sites());
    for (unsigned i=0; i<all_sites.size(); ++i) all_sites[i] = i;
    std::shuffle(all_sites.begin(), all_sites.end(), Simulator::rng);
    // first 'num_doped_sites' of these shuffled sites are doped
    // bonds connecting to these sites are impurity bonds
    out_bond_iterator ob, ob_end;
    in_bond_iterator ib, ib_end;
    unsigned impurity_btype = Simulator::Model::get_impurity_bondtype();
    for (unsigned i=0; i<num_doped_sites; ++i) {
      for (std::tie(ob,ob_end)=out_bonds(all_sites[i]); ob!=ob_end; ++ob) {
        Simulator::LatticeGraph::change_type_value(ob, impurity_btype);
      }
      for (std::tie(ib,ib_end)=in_bonds(all_sites[i]); ib!=ib_end; ++ib) {
        Simulator::LatticeGraph::change_type_value(ib, impurity_btype);
      }
    }
  }
} 

int MCSimulation::start(input::Parameters& parms)
{
  // update model parameters
  Simulator::update_parameters(parms);

  // observables
  observables.reset();

  /*-----------------Simulation START-----------------*/
  // initial state
  Simulator::init_state_random();
  // warmup runs
  for (int i=0; i<warmup; ++i) Simulator::update_state_metropolis();
  // measurement
  int count = 0;
  int num_measurement = 0;
  while (num_measurement != measure_samples) {
    Simulator::update_state_metropolis();
    if (count++ == min_interval) {
      // make measurements
      count = 0;
      num_measurement++;
      do_measurements();
    }
  }
  // output
  int z = 6;
  double J = parms.set_value("J", 1.0);
  double T_star = Simulator::kB * Simulator::T/(z*J);
  observables.print(T_star);
  /*-----------------Simulation END-----------------*/
  return 0;
}

// measurements
inline void MCSimulation::do_measurements(void)
{
  // energy
  if (need_energy) {
    energy_terms = Simulator::get_energy();
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
    double x = Simulator::get_magnetization();
    if (observables.magn()) observables.magn() << x;
    if (observables.magn_sq()) observables.magn_sq() << x*x;
  }

  // Potts magnetization
  if (observables.potts_magn() || observables.potts_magn_sq()) {
    double x = Simulator::get_potts_magnetization();
    if (observables.potts_magn()) observables.potts_magn() << x;
    if (observables.potts_magn_sq()) observables.potts_magn_sq() << x*x;
  }

  // strain
  if (observables.strain() || observables.strain_sq()) {
    double x = Simulator::get_strain();
    if (observables.strain()) observables.strain() << x;
    if (observables.strain_sq()) observables.strain_sq() << x*x;
  }
}


