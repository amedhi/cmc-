/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MCSIMULATION_H
#define MCSIMULATION_H
#include <map>
#include <string>
#include <cmc/simulation.h>

class MCSimulation : public mc::Simulator
{
public:
  MCSimulation(input::Parameters& parms);
  ~MCSimulation() {}
  int start(input::Parameters& parms);
private:
  void do_measurements(void);
  std::map<std::string, double> print_parms;
};



#endif
