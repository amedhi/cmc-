/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include <iostream>
#include <string>
#include <cmc/simulation.h>

/*------Your Simulation class to be derived from mc::Simulator-----*/
class MCSimulation : public mc::Simulator
{
public:
  MCSimulation(input::Parameters& parms);
  ~MCSimulation() {}
  int start(input::Parameters& parms);
private:
  void do_measurements(void);
};

/*-----Implementation----*/
// Constructor
MCSimulation::MCSimulation(input::Parameters& parms) : Simulator(parms) 
{
  // You need observables as function of what?
  observables.as_function_of("T");

  // Define your observable operators (other than energy)
  if (observables.magn() || observables.magn_sq() ) {
    Simulator::set_magn_op("S(i)");
  }

  /* 
   * If you need the value of some input parameter (say, "rng_seed") 
   * from the input file, you can get it like this:
   */
   int seed = parms.set_value("rng_seed", 0); 
   // 0 is the default value set in case the parameter "rng_seed" is not found
   // in the input file.
   // std::cout << "rng_seed = " << seed << std::endl;
} 

// Actual simulation inside this function
int MCSimulation::start(input::Parameters& parms)
{
  // update simulation parameters
  Simulator::update_parameters(parms);

  // observable values must be reset before simulation start
  observables.reset();

  /*-----------------Simulation START-----------------*/
  // initial state
  Simulator::init_state_random();
  // warmup runs
  for (int i=0; i<warmup; ++i) Simulator::update_state();
  // production runs
  int count = 0;
  int num_measurement = 0;
  while (num_measurement != measure_samples) {
    Simulator::update_state();
    if (count++ == min_interval) {
      // make measurements & store the results
      count = 0;
      num_measurement++;
      do_measurements();
    }
  }
  // print the results for calculated observables
  observables.print(T);
  /*-----------------Simulation END-----------------*/
  // done!
  return 0;
}

// Measurement of observables
inline void MCSimulation::do_measurements(void)
{
  // energy
  if (Simulator::need_energy) {
    energy_terms = Simulator::get_energy();
    // total energy
    if (observables.energy()) {
      double e = energy_terms.sum();
      observables.energy() << e;
    }
    // [total energy]^2
    if (observables.energy_sq()) {
      double e = energy_terms.sum();
      observables.energy_sq() << e*e;
    }
    // energy for individual Hamiltonian terms
    if (observables.energy_terms()) {
      observables.energy_terms() << energy_terms;
    }
    // [energy for individual Hamiltonian terms]^2
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
}

// In 'main' just feed your Simulatior to the scheduler
int main(int argc, const char *argv[])
{
  try {
   return scheduler::start(argc, argv, scheduler::Task<MCSimulation>());
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }
  catch (...) {
    std::cout << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
}



