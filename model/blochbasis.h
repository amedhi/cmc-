/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-13 23:41:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-18 01:10:25
*----------------------------------------------------------------------------*/
#ifndef BLOCHBASIS_H
#define BLOCHBASIS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <stdexcept>
#include "../lattice/lattice.h"
#include "../lattice/graph.h"

namespace basis {

using Vector3d = Eigen::Vector3d; 
using Vector3i = Eigen::Vector3i; 
using basis_state = lattice::graph::vertex_descriptor;

class BlochBasis 
{
public:
  // ctors
  BlochBasis() {}
  BlochBasis(const lattice::Lattice& lattice, const lattice::graph::LatticeGraph& graph)
    { construct(lattice, graph); }
  ~BlochBasis() {}
  // setter functions
  int construct(const lattice::Lattice& lattice, const lattice::graph::LatticeGraph& graph);
  // getter functions
  unsigned dimension(void) const { return basis_states.size(); }
  unsigned num_kpoints(void) const { return kpoints.size(); }
  Vector3d bloch_vector(const unsigned& k) const { return kpoints[k]; }
  basis_state representative_state(const basis_state& s, const lattice::graph::LatticeGraph& graph,
    Vector3d& R) const;
  basis_state site_basis(const unsigned& idx) const 
    { return basis_states[idx]; }
  unsigned state_idx(const basis_state& s) const 
    {  auto pos = state_indices.find(s); return pos != state_indices.end() ? pos->second : null_index; }
  unsigned null_idx(void) const { return null_index; }
  //double MatrixElementent(const op_iterator& op_it, const unsigned& type) { return op_it->second[type].value; }
  // friends
private:
  Vector3d b1, b2, b3;
  std::vector<Vector3d> kpoints;
  //std::vector<Vector3i> translation_vectors;
  std::vector<basis_state> basis_states;
  //std::vector<unsigned> basis_state;
  std::unordered_map<basis_state, unsigned> state_indices;
  unsigned null_index;

  // helper functions
  void make_bloch_vectors(const lattice::Lattice& lattice);
  void make_site_basis(const lattice::Lattice& lattice, const lattice::graph::LatticeGraph& graph);
};


} // end namespace basis

#endif
