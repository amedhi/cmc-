/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: amedhi
* Date:   2016-03-01 00:11:01
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-17 23:52:15
*----------------------------------------------------------------------------*/
#include "graph.h"

namespace lattice {
namespace graph {

//using namespace boost;
//LatticeGraph::LatticeGraph(const Lattice& lattice) 
//{
//  construct(lattice);
//}

LatticeGraph::LatticeGraph(input::Parameters& parms) : lattice_(parms)
{
  construct_graph();
} 

void LatticeGraph::construct(input::Parameters& parms) 
{
  lattice_.construct(parms);
  construct_graph();
} 

void LatticeGraph::construct_graph(void) 
{
  // all the sites and the bonds
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (unsigned i=0; i<lattice_.num_unitcells(); ++i) {
    translated_cell = lattice_.get_translated_cell(bravindex);
    // collect the sites
    for (unsigned n=0; n<translated_cell.num_sites(); ++n) sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (unsigned n=0; n<translated_cell.num_bonds(); ++n) {
      Bond b = translated_cell.bond(n);
      if (lattice_.connect_bond(b)) {
        //b.reset_id(bonds.size());
        bonds.push_back(b);
      }
    }
    bravindex = lattice_.get_next_bravindex(bravindex);
  }

  // construct the graph
  clear();
  vertex_types_set_.clear();

  // add vertices 
  for (unsigned i=0; i<sites.size(); ++i) {
    vertex_descriptor v = add_vertex(*this);
    // vertex properties
    this->operator[](v).uid = sites[i].uid();
    this->operator[](v).type = sites[i].type();
    this->operator[](v).stype = sites[i].type();
    this->operator[](v).atomid = sites[i].atomid();
    this->operator[](v).bravindex = sites[i].bravindex();
    this->operator[](v).coord = sites[i].coord();
    this->operator[](v).cell_coord = sites[i].cell_coord();

    // collect type value in a set
    vertex_types_set_.insert(sites[i].type());
  } 
  num_vertices_ = boost::num_vertices(*this);
  boost::tie(vi_begin_, vi_end_) = boost::vertices(*this);
  vertex_index_map = boost::get(boost::vertex_index, *this);
  // is the vertex_descriptor same as its 'position' or 'id'?
  for (unsigned i=0; i<num_vertices_; ++i) {
    if (i != boost::vertex(i, *this)) {
      throw std::logic_error("LatticeGraph::construct_graph: vertex_descriptor assumption wrong");
    }
  } 

  //vertex_uid_map = get(&VertexProperties::uid, g); 
  //vertex_cellcord_map = get(&VertexProperties::cell_coord, g); 

  // add the edges
  edge_descriptor e;
  bool flag;
  vertex_descriptor source, target;
  for (unsigned i=0; i<bonds.size(); ++i) {
    source = boost::vertex(bonds[i].src_id(), *this);
    target = boost::vertex(bonds[i].tgt_id(), *this);
    boost::tie(e, flag) = add_edge(source, target, *this);
    //std::cout << e << "\n";
    // edge properties
    this->operator[](e).type = bonds[i].type();
    this->operator[](e).stype = bonds[i].type();
    this->operator[](e).sign = bonds[i].sign();
    this->operator[](e).vector = sites[bonds[i].tgt_id()].coord() - sites[bonds[i].src_id()].coord();
    //std::cout << "bond: " << vertex_index_map[source] << " --- " << vertex_index_map[target] << "\n";
    //std::cout << bonds[i].src_offset() << "\n";
    //std::cout << bonds[i].tgt_offset() << "\n\n";
  }
  num_edges_ = boost::num_edges(*this);
  boost::tie(ei_begin_, ei_end_) = boost::edges(*this);
  //for (unsigned i=0; i<num_edges_; ++i) {
  //  std::cout << boost::edge(i, *this) << "\n"; 
  //}
  //std::cout << "Vertices = " << num_vertices(g) << std::endl;
  //std::cout << "Edges = " << num_edges(g) << std::endl;
  //std::pair<adjacency_iterator,adjacency_iterator> p = boost::adjacent_vertices(0, *this);
  //for (auto it=p.first; it!=p.second;++it) { std::cout << *it << "\n"; }
  //std::pair<in_edge_iterator, in_edge_iterator> po = boost::in_edges(0, *this); 
  //for (auto it=po.first; it!=po.second;++it) { std::cout << *it << "\n"; }
}

} // end namespace graph
} // end namespace lattice