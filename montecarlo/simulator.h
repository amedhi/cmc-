/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <Eigen/Core>
#include "../scheduler/worker.h"
#include "./random.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "sitebasisstate.h"
#include "observables.h"
#include "observable_operator.h"

namespace mc {

//class Simulator : public std::vector<SiteState>
class Simulator : public lattice::graph::LatticeGraph, public model::Model, 
                  public scheduler::Worker 
{
public:
  using SiteState = SiteBasisState;
  using state_idx = SiteBasisState::state_idx;
  //Simulator() {};
  Simulator(input::Parameters& parms); 
  ~Simulator() {}
  void start(input::Parameters& parms);
  void run(void) {} 
  void finish(void) {} 
  void dostep(void) {} 
  void halt(void) {} 
  static void print_copyright(std::ostream& os);
private:
  RandomNumber rng;
  std::vector<SiteBasisState> state;
  std::array<Eigen::MatrixXd, lattice::MAX_BOND_TYPES> boltzmann_table;

  // system parameters
  double T; // temperature
  double kB; 
  double beta; 

  // mc parameters
  int measure_samples; 
  int warmup;
  int min_interval;

  // observables
  Observables observables;

  void init(void);
  void init_state_random(void);
  void init_boltzmann_table(void);
  void update_state_metropolis(void);
  void do_measurements(void);
  mc::VectorData get_energy(void);
  double get_magnetization(void);
  double get_potts_magnetization(void);
  double get_strain(void);
};


} // end namespace monte carlo

#endif
