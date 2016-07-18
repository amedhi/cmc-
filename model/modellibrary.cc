/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-17 23:51:58
*----------------------------------------------------------------------------*/
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Model::define_model(const lattice::Lattice& lattice, const input::Parameters& inputs)
{
  //int info;
  unsigned ntypes;
  std::vector<MatrixElement> matrix_elem(20);

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  /*------------- 'Hubbard' model--------------*/
  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    switch (lattice.id()) {
      /*------------- 'SQUARE' lattice--------------*/
      case lattice::lattice_id::SQUARE:
        // model parameters
        set_parameter(param::t, "t", inputs);
        set_parameter(param::U, "U", inputs);

        // hopping operators
        matrix_elem[0] = product(-1.0, param::t);
        add_operator(op_id::upspin_hop, matrix_elem, ntypes=1);
        add_operator(op_id::dnspin_hop, matrix_elem, ntypes=1);

        break;
      default:
        throw std::range_error("*error: modellibrary: model not defined for the given lattice");
    }
  }

  /*------------- 't-J' model--------------*/
  else if (model_name == "TJ") {
    mid = model_id::TJ;
    switch (lattice.id()) {
      /*------------- 'SQUARE' lattice--------------*/
      case lattice::lattice_id::SQUARE:
        // model parameters
        set_parameter(param::t, "t", inputs);
        set_parameter(param::J, "J", inputs);

        // hopping operators
        matrix_elem[0] = product(-1.0, param::t);
        matrix_elem[2] = product(-1.0, param::t);
        add_operator(op_id::upspin_hop, matrix_elem, ntypes=2);
        add_operator(op_id::dnspin_hop, matrix_elem, ntypes=2);

        break;
      default: 
        throw std::range_error("*error: modellibrary: model not defined for the given lattice");
    }
  }

  /*------------- undefined lattice--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  return 0;
}

int Model::define_mf_model(const lattice::Lattice& lattice, const input::Parameters& inputs)
{
  // Mean-field version

  //int info;
  unsigned ntypes;
  std::vector<MatrixElement> matrix_elem(20);

  // input information
  read_mf_order(inputs);

  // first copy the 'single particle operators'
  for (auto it=operators.begin(); it!=operators.end(); ++it) {
    if (it->second.single_particle_op) mf_operators.insert({*it});
  }

  // next the mean-field order terms
  // mean-field order
  switch (mf_order) {
    case mforder::SC:
      set_parameter(param::delta_sc, "delta_sc", inputs);
      switch (lattice.id()) {
        case lattice::lattice_id::SQUARE:
          switch (pairing_symmetry) {
            case pairsymm::extended_s:
              matrix_elem[0] = product(+1.0, param::delta_sc);
              matrix_elem[1] = product(+1.0, param::delta_sc);
              break;
            case pairsymm::dwave:
              matrix_elem[0] = product(+1.0, param::delta_sc);
              matrix_elem[1] = product(-1.0, param::delta_sc);
              break;
            default: break;
          }
          add_operator(op_id::bondsinglet_hop, matrix_elem, ntypes=2);
          break;
        default: break;
      }
      break;
    default: break;
  }
  return 0;
}

int Model::read_mf_order(const input::Parameters& inputs)
{
  int info;
  mf_order_name = inputs.set_value("mf_order", "NONE", info);
  boost::to_upper(mf_order_name);
  if (mf_order_name == "NONE") {
    mf_order_name = "FERMISEA";
    mf_order = mforder::Fermisea;
  }
  else if (mf_order_name == "FERMISEA") {
    mf_order = mforder::Fermisea;
  }
  else if (mf_order_name == "AF") {
    mf_order = mforder::AF;
  }
  else if (mf_order_name == "FM") {
    mf_order = mforder::FM;
  }
  else if (mf_order_name == "SC") {
    mf_order = mforder::SC;
    read_pairing_symmetry(inputs);
  }
  else if (mf_order_name == "TRIPLET_SC") {
    mf_order = mforder::Triplet_SC;
    read_pairing_symmetry(inputs);
  }
  else {
    throw std::range_error("*error: modellibrary: undefined 'mf_order'");
  }
  return 0;
}

int Model::read_pairing_symmetry(const input::Parameters& inputs)
{
  int info;
  pairing_symmetry_name = inputs.set_value("pairing_symmetry", "NONE", info);
  boost::to_upper(pairing_symmetry_name);
  if (pairing_symmetry_name == "SWAVE") {
    pairing_symmetry = pairsymm::swave;
  }
  else if (pairing_symmetry_name == "EXTENDED_S") {
    pairing_symmetry = pairsymm::extended_s;
  }
  else if (pairing_symmetry_name == "DWAVE") {
    pairing_symmetry = pairsymm::dwave;
  }
  else if (pairing_symmetry_name == "D_PLUS_ID") {
    pairing_symmetry = pairsymm::d_plus_id;
  }
  else if (pairing_symmetry_name == "PWAVE") {
    pairing_symmetry = pairsymm::pwave;
  }
  else if (pairing_symmetry_name == "P_PLUS_IP") {
    pairing_symmetry = pairsymm::p_plus_ip;
  }
  else {
    throw std::range_error("*error: modellibrary: undefined 'pairing_symmetry'");
  }
  return 0;
}

int Model::construct(const lattice::Lattice& lattice, const input::Parameters& inputs)
{
  // reset
  parameters.clear();
  operators.clear();
  mf_operators.clear();

  define_model(lattice, inputs);
  define_mf_model(lattice, inputs);
  finalize();

  return 0;
}


} // end namespace model
