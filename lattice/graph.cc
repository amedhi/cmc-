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

LatticeGraph::LatticeGraph(const Lattice& lattice) 
{
  construct(lattice);
}

void LatticeGraph::construct(const Lattice& lattice) 
{
  // all the sites and the bonds
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (unsigned i=0; i<lattice.num_unitcells(); ++i) {
    translated_cell = lattice.get_translated_cell(bravindex);
    // collect the sites
    for (unsigned n=0; n<translated_cell.num_sites(); ++n) sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (unsigned n=0; n<translated_cell.num_bonds(); ++n) {
      Bond b = translated_cell.bond(n);
      if (lattice.connect_bond(b)) {
        //b.reset_id(bonds.size());
        bonds.push_back(b);
      }
    }
    bravindex = lattice.get_next_bravindex(bravindex);
  }

  // construct the graph
  g.clear();

  // add vertices 
  for (unsigned i=0; i<sites.size(); ++i) {
    graph::vertex_descriptor v = add_vertex(g);
    // vertex properties
    g[v].uid = sites[i].uid();
    g[v].type = sites[i].type();
    g[v].stype = sites[i].type();
    g[v].atomid = sites[i].atomid();
    g[v].bravindex = sites[i].bravindex();
    g[v].coord = sites[i].coord();
    g[v].cell_coord = sites[i].cell_coord();
  } 
  //boost::tie(vi_begin, vi_end) = vertices(g);
  vertex_index_map = get(boost::vertex_index, g);
  //vertex_uid_map = get(&VertexProperties::uid, g); 
  //vertex_cellcord_map = get(&VertexProperties::cell_coord, g); 

  // add the edges
  graph::edge e;
  bool flag;
  graph::vertex_descriptor source, target;
  for (unsigned i=0; i<bonds.size(); ++i) {
    source = boost::vertex(bonds[i].src_id(), g);
    target = boost::vertex(bonds[i].tgt_id(), g);
    boost::tie(e, flag) = add_edge(source, target, g);
    // edge properties
    g[e].type = bonds[i].type();
    g[e].stype = bonds[i].type();
    g[e].sign = bonds[i].sign();
    g[e].vector = sites[bonds[i].tgt_id()].coord() - sites[bonds[i].src_id()].coord();
    //std::cout << "bond: " << vertex_index[source] << " --- " << vertex_index[target] << "\n";
    //std::cout << bonds[i].src_offset() << "\n";
    //std::cout << bonds[i].tgt_offset() << "\n\n";
  }
  boost::tie(ei_begin, ei_end) = edges(g);
  //std::cout << "Vertices = " << num_vertices(g) << std::endl;
  //std::cout << "Edges = " << num_edges(g) << std::endl;
}

} // end namespace graph
} // end namespace lattice
