/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-15 14:51:39
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "../lattice/lattice.h"
#include "model.h"
#include "blochbasis.h"

namespace model {

using Matrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
using Complex = std::complex<double>;

class Hamiltonian 
{
public:
  // ctors
  //Hamiltonian() {}
  Hamiltonian(const lattice::Lattice& lattice, const lattice::graph::LatticeGraph& graph, 
    const model::Model& model);
  ~Hamiltonian() {}
  // setter functions
  unsigned num_blocks(void) { return basis.num_kpoints(); }
  int construct_block(const unsigned& k, const model::Model& model, const lattice::graph::LatticeGraph& graph);
  // getter functions
  //double MatrixElementent(const op_iterator& op_it, const unsigned& type) { return op_it->second[type].value; }
  // friends
  //double coupling_constant(const op_iterator& op, const unsigned& type) { return upspin_hop[type]; }
private:
   basis::BlochBasis basis;
   unsigned dim;
   Matrix h0_up, h0_down;
   Matrix hk_up, hk_down;
   Matrix bdg;
  //Eigen::Matrix hk;
  //Eigen::Matrix bdg;

   // helper functions
   int construct_h0(const model::Model& model, const lattice::graph::LatticeGraph& graph);
   int construct_hk(const unsigned& k, const model::Model& model, const lattice::graph::LatticeGraph& graph);
};


} // end namespace model

#endif