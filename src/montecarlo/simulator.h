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
#include "../observable/observables.h"
#include "../observable/observable_operator.h"

namespace mc {

class Simulator : public lattice::graph::LatticeGraph, public model::Model, 
                  public scheduler::Worker 
{
public:
  using SiteState = SiteBasisState;
  using state_idx = SiteBasisState::state_idx;
  //Simulator() {};
  Simulator(input::Parameters& parms); 
  ~Simulator() {}
  using Model::update_parameters;
  void start(input::Parameters& parms) override;
  void run(void) override {} 
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);

protected:
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
  virtual void do_measurements(void);
  virtual mc::VectorData get_energy(void);
  virtual double get_magnetization(void);
  virtual double get_potts_magnetization(void);
  virtual double get_strain(void);
};


} // end namespace monte carlo

#endif
