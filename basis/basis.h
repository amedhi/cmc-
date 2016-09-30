/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <string>
#include <vector>
#include "qn.h"

namespace basis {

class State
{
  private:
    QN_set qn_;
    std::size_t idx_;
  public:
    State(): idx_(0) {}
    State(size_t idx, QN_set qn): qn_(qn), idx_(idx) {}
    State(const State &s): qn_(s.qn_),idx_(s.idx_) {}
};

class Basis : public std::vector<State>
{
  private:
    short qn;
public:
  Basis();
  ~Basis();
};



} // end namespace basis

#endif
