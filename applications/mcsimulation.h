/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MCSIMULATION_H
#define MCSIMULATION_H
#include <montecarlo/simulator.h>
//#include <cmc++/cmc.h>

class MCSimulation : public mc::Simulator
{
public:
  MCSimulation(input::Parameters& parms);
  ~MCSimulation() {}
  void start(input::Parameters& parms);
private:
  void do_measurements(void);
};



#endif
