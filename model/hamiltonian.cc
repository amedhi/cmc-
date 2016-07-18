/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-15 14:51:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-17 23:51:52
*----------------------------------------------------------------------------*/
#include "hamiltonian.h"

namespace model {

Hamiltonian::Hamiltonian(const lattice::Lattice& lattice, const lattice::graph::LatticeGraph& graph, 
    const model::Model& model) : basis(lattice, graph)
{
  // initialization
  dim = basis.dimension();
  h0_up.resize(dim, dim);
  hk_up.resize(dim, dim);
  h0_down.resize(dim, dim);
  hk_down.resize(dim, dim);
  bdg.resize(dim, dim);

  // hamiltoian matrix for intracell terms 
  construct_h0(model, graph);
}

int Hamiltonian::construct_block(const unsigned& k, const model::Model& model, 
  const lattice::graph::LatticeGraph& graph)
{
  construct_hk(k, model, graph);
  return 0;
}

int Hamiltonian::construct_hk(const unsigned& k, const model::Model& model, 
  const lattice::graph::LatticeGraph& graph)
{
  model::op_iterator op, op_end;
  lattice::graph::out_edge_iterator ei, ei_end;
  lattice::graph::vertex_descriptor s, t, rs;
  Complex term;
  basis::Vector3d kvec, R;
  unsigned i, j, type;

  // k-vector
  kvec = basis.bloch_vector(k);

  //  spin-UP block 
  hk_up = h0_up;
  // upspin_hop
  for (std::tie(op, op_end)=model.mf_operators_upspin_hop(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      //std::cout << "s = " << s  << "\n";
      for (std::tie(ei, ei_end)=graph.out_edges(s); ei != ei_end; ++ei) {
        //std::cout << "edge type = " << graph.edge_type(ei) << "\n";
        t = graph.target_vertex(ei);
        j = basis.state_idx(t);
        if (j == basis.null_idx()) {
          type = graph.edge_type(ei);
          term = model.matrix_element(op, type);
          rs = basis.representative_state(t, graph, R);
          term = term * std::exp(+ii() * kvec.dot(R));
          j = basis.state_idx(rs);
          hk_up(i,j) += term;
          hk_up(j,i) += std::conj(term);
        }
      }
    }
  }
  //std::cout << hk_up << "\n";

  //  spin-DOWN block 
  if (model.has_spinflip_symmetry()) return 0; 
  hk_down = h0_down;
  // upspin_hop
  for (std::tie(op, op_end)=model.mf_operators_dnspin_hop(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      //std::cout << "s = " << s  << "\n";
      for (std::tie(ei, ei_end)=graph.out_edges(s); ei != ei_end; ++ei) {
        //std::cout << "edge type = " << graph.edge_type(ei) << "\n";
        t = graph.target_vertex(ei);
        j = basis.state_idx(t);
        if (j == basis.null_idx()) {
          type = graph.edge_type(ei);
          term = model.matrix_element(op, type);
          rs = basis.representative_state(s, graph, R);
          term = term * std::exp(+ii() * kvec.dot(R));
          j = basis.state_idx(rs);
          hk_down(i,j) += term;
          hk_down(j,i) += std::conj(term);
        }
      }
    }
  }
  //std::cout << hk_down << "\n";

  return 0;
}

int Hamiltonian::construct_h0(const model::Model& model, const lattice::graph::LatticeGraph& graph) 
{
  model::op_iterator op, op_end;
  lattice::graph::out_edge_iterator ei, ei_end;
  lattice::graph::vertex_descriptor s, t;
  Complex term;
  unsigned i, j, type;

  //  spin-UP block 
  h0_up.setZero();
  // upspin_hop
  for (std::tie(op, op_end)=model.mf_operators_upspin_hop(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      //std::cout << "s = " << s  << "\n";
      for (std::tie(ei, ei_end)=graph.out_edges(s); ei != ei_end; ++ei) {
        //std::cout << "edge type = " << graph.edge_type(ei) << "\n";
        t = graph.target_vertex(ei);
        j = basis.state_idx(t);
        if (j != basis.null_idx()) {
          type = graph.edge_type(ei);
          term = model.matrix_element(op, type);
          h0_up(i,j) += term;
          h0_up(j,i) += std::conj(term);
        }
      }
    }
  }
  // ni_up
  for (std::tie(op, op_end)=model.mf_operators_niup(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      type = graph.vertex_type(s);
      term = model.matrix_element(op, type);
      h0_up(i,i) += term;
    }
  }
  //std::cout << h0_up << "\n";

  //  spin-DOWN block 
  if (model.has_spinflip_symmetry()) return 0; 
  h0_down.setZero();
  // dnspin_hop
  for (std::tie(op, op_end)=model.mf_operators_dnspin_hop(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      //std::cout << "s = " << s  << "\n";
      for (std::tie(ei, ei_end)=graph.out_edges(s); ei != ei_end; ++ei) {
        //std::cout << "edge type = " << graph.edge_type(ei) << "\n";
        t = graph.target_vertex(ei);
        j = basis.state_idx(t);
        if (j != basis.null_idx()) {
          type = graph.edge_type(ei);
          term = model.matrix_element(op, type);
          h0_down(i,j) += term;
          h0_down(j,i) += std::conj(term);
        }
      }
    }
  }
  // ni_down
  for (std::tie(op, op_end)=model.mf_operators_nidn(); op!=op_end; ++op) {
    for (i=0; i<dim; ++i) {
      s = basis.site_basis(i);
      type = graph.vertex_type(s);
      term = model.matrix_element(op, type);
      h0_down(i,i) += term;
    }
  }
  //std::cout << h0_down << "\n";

  return 0;
}

//void Hamiltonian::h0(const lattice::graph::LatticeGraph& graph, const model::Model& model) 
//{}


} // end namespace model
