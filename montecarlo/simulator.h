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
  // typedefs
  //using site_iterator = LatticeGraph::site_iterator;
  //using bond_iterator = LatticeGraph::bond_iterator;
  //using site_descriptor = LatticeGraph::site_descriptor;
  //using siteterm_iterator = Model::siteterm_iterator;
  //using bondterm_iterator = Model::bondterm_iterator;
  // system state
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


  //alps::RealObservable Magnetization("Magnetization");

  void init(void);
  void init_state_random(void);
  void init_boltzmann_table(void);
  void update_state_metropolis(void);
  void init_observables(const input::Parameters& parms);
  void do_measurements(void);
  mc::VectorData get_energy(void);
  double magnetization(void);
  double potts_magnetization(void);
  double strain(void);
};


} // end namespace monte carlo

#endif
