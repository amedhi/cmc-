/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: amedhi
* Date:   2016-03-01 00:11:01
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-18 01:08:36
*----------------------------------------------------------------------------*/
#ifndef GRAPH_H
#define GRAPH_H

#include "lattice.h"
#include <boost/graph/adjacency_list.hpp>
//#include <Eigen/Dense>

namespace lattice {
namespace graph {

// Vertex properties
struct VertexProperties {
  unsigned uid; // unitcell id
  unsigned type; 
  unsigned stype; // symmetry type
  unsigned atomid; 
  Vector3i bravindex;
  Vector3d coord; 
  Vector3d cell_coord;
};

// Edge properties
struct EdgeProperties {
  unsigned type; 
  unsigned stype; 
  int sign; 
  Vector3d vector; 
};

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, 
              VertexProperties, EdgeProperties>;
using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
using edge = boost::graph_traits<Graph>::edge_descriptor;
using vertex_iterator = boost::graph_traits<Graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<Graph>::edge_iterator;
using out_edge_iterator = boost::graph_traits<Graph>::out_edge_iterator;


class LatticeGraph 
{
public:
  // ctors
  LatticeGraph() { g.clear(); }
  LatticeGraph(const Lattice& lattice);
  ~LatticeGraph() {}

  // setter functions
  void construct(const Lattice& lattice);

  // getter functions
  unsigned num_vertices(void) const { return boost::num_vertices(g); }
  vertex_iterator vertex_begin(void) const { return vi_begin; }
  vertex_iterator vertex_end(void) const { return vi_end; }
  //unsigned vertex_id(const vertex_iterator& vi) const { return vertex_index_map[*vi]; }
  //unsigned vertex_id(const unsigned& i) const { return vertex_index_map[boost::vertex(i, g)]; }
  std::pair<vertex_iterator, vertex_iterator> vertices(void) const { return boost::vertices(g); }
  vertex_descriptor vertex(const unsigned& i) const { return boost::vertex(i, g); }
  unsigned vertex_id(const vertex_descriptor& v) const { return vertex_index_map[v]; }
  unsigned vertex_uid(const vertex_descriptor& v) const { return g[v].uid; }
  Vector3d vertex_cellcord(const vertex_descriptor& v) const { return g[v].cell_coord; }
  //unsigned vertex_uid(const unsigned& i) const { return g[boost::vertex(i, g)].uid; }
  unsigned vertex_type(const vertex_descriptor& v) const { return g[v].type; }
  unsigned vertex_type(const vertex_iterator& vi) const { return g[*vi].type; }
  unsigned edge_type(const out_edge_iterator& ei) const { return g[*ei].type; }
  vertex_descriptor target_vertex(const out_edge_iterator& ei) const { return target(*ei, g); }
  //using out_edges = boost::out_edges;
  std::pair<out_edge_iterator, out_edge_iterator> out_edges(const vertex_descriptor& v) const 
    { return boost::out_edges(v, g); }
  // friends
private:
  Graph g;
  vertex_iterator vi_begin, vi_end;
  edge_iterator ei_begin, ei_end;
  boost::property_map<graph::Graph, boost::vertex_index_t>::type vertex_index_map;
  //boost::property_map<graph::Graph, unsigned VertexProperties::*>::type vertex_uid_map;
  //boost::property_map<graph::Graph, Vector3d VertexProperties::*>::type vertex_cellcord_map;
};


} // end namespace graph
} // end namespace lattice

#endif
