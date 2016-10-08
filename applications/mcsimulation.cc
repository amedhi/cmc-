/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "mcsimulation.h"

MCSimulation::MCSimulation(input::Parameters& parms) : Simulator(parms) 
{
} 

void MCSimulation::start(input::Parameters& parms)
{
  // update model parameters
  Simulator::update_parameters(parms);

  // set variable parameters
  T = parms.set_value("T", 1.0);
  kB = parms.set_value("kB", 1.0);
  beta = 1.0/(kB*T);

  // values of exp(-beta*E) for bond energy (E) for Metropolis algorithm
  Simulator::init_boltzmann_table();

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
  observables.print(T);
  /*-----------------Simulation END-----------------*/
}

// measurements
inline void MCSimulation::do_measurements(void)
{
  // energy
  if (observables.energy_terms()) {
    observables.energy_terms() << get_energy();
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
    if (observables.strain_sq()) observables.strain_sq() << x*x;
  }
}


