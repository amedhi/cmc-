/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 12:08:30
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-18 00:47:24
*----------------------------------------------------------------------------*/
#include <iostream>
#include "scheduler/task.h"
#include "lattice/lattice.h"
#include "lattice/graph.h"
#include "model/model.h"
#include "model/blochbasis.h"
#include "model/hamiltonian.h"
#include "mytask.h"

void MyTask::start(input::Parameters& parms)
{
  //double u = parms.set_value("U", 1.0);
  std::cout << " task " << parms.task_id()+1 << " of " << parms.task_size() << std::endl;

  std::cout << "\nChecking latticelibrary" << std::endl;
  lattice::Lattice lattice(parms); 
  lattice::graph::LatticeGraph graph(lattice);
  //basis::BlochBasis bloch_basis(lattice, graph);
  model::Model model(lattice, parms);
  model::Hamiltonian ham(lattice, graph, model);
  for (unsigned k=0; k<ham.num_blocks(); ++k) {
    ham.construct_block(k, model, graph);
  }
  //hamiltonian.make_block(block k, eigen_value)

  /*lattice::graph::vertex_iterator vi;
  lattice::graph::edge_iterator ei;
  for (vi=graph.vertex_begin(); vi!=graph.vertex_end(); ++vi) {
    std::cout << "id = " << graph.vertex_id(vi) << std::endl;
    std::cout << "type = " << graph.vertex_type(vi) << std::endl;
  }*/
}
