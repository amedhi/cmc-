/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:46
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <stdexcept>
#include "../lattice/lattice.h"

namespace model {

/*---------------model names-----------------*/
enum class model_id {
  UNDEFINED, HUBBARD, TJ
};

// parameter names
enum class param {
  t, t2, t3, th, e0, J, J2, J3, Jh, U,
  delta_af, delta_fm, delta_sc, delta_00, delta_10, delta_11, delta_m1
};

enum class spin {UP, DOWN, UD};

// operators
enum class op_id {
  upspin_hop, dnspin_hop, ni_up, ni_dn, ni, Hubbard, Exchange, bondsinglet_hop
};

// mean-field order
enum class mforder {
  Fermisea, AF, FM, SC, Triplet_SC 
};

// SC pairing symmetry
enum class pairsymm {
  undefined, swave, extended_s, dwave, d_plus_id, pwave, p_plus_ip 
};

using Complex = std::complex<double>;
struct MatrixElement {Complex value; Complex factor; param pname;};
struct op_value {spin sigma; std::vector<MatrixElement> matrix_elem; bool single_particle_op;};

// iterators
using op_iterator = std::multimap<op_id, op_value>::const_iterator;

class Model 
{
public:
  // ctors
  Model() {}
  Model(const lattice::Lattice& lattice, const input::Parameters& inputs) { construct(lattice, inputs); }
  ~Model() {}
  // setter functions
  int construct(const lattice::Lattice& lattice, const input::Parameters& inputs);
  int define_model(const lattice::Lattice& lattice, const input::Parameters& inputs);
  int define_mf_model(const lattice::Lattice& lattice, const input::Parameters& inputs);
  // getter functions
  std::pair<op_iterator, op_iterator> mf_operators_upspin_hop(void) const 
    { return {mf_ciup_begin, mf_ciup_end}; }
  std::pair<op_iterator, op_iterator> mf_operators_dnspin_hop(void) const 
    { return {mf_cidn_begin, mf_cidn_end}; }
  std::pair<op_iterator, op_iterator> mf_operators_niup(void) const 
    { return {mf_niup_begin, mf_niup_end}; }
  std::pair<op_iterator, op_iterator> mf_operators_nidn(void) const 
    { return {mf_nidn_begin, mf_nidn_end}; }

  op_iterator mf_operator_begin(void) const { return mf_operators.begin(); }
  op_iterator mf_operator_end(void) const { return mf_operators.end(); }
  op_id operator_id(const op_iterator& it) const { return it->first; }
  Complex matrix_element(const op_iterator& it, const unsigned& type) const 
    { return it->second.matrix_elem[type].value; }
  bool has_spinflip_symmetry(void) const { return spinflip_symmetry; }
  bool has_inversion_symmetry(void) const { return inversion_symmetry; }
  //double MatrixElementent(const op_iterator& op_it, const unsigned& type) { return op_it->second[type].value; }
  // friends
  //double coupling_constant(const op_iterator& op, const unsigned& type) { return upspin_hop[type]; }
private:

  // model-id
  model_id mid {model_id::UNDEFINED};
  std::string model_name {""};
  //
  bool spinflip_symmetry {true};
  bool inversion_symmetry {true};
  // parameters
  std::map<param, double> parameters; 
  // operators
  std::multimap<op_id, op_value> operators;

  // mean-field model
  std::string mf_order_name;
  mforder mf_order = mforder::Fermisea;
  pairsymm pairing_symmetry = pairsymm::undefined;
  std::string pairing_symmetry_name;
  // mean-field operators
  std::multimap<op_id, op_value> mf_operators;
  op_iterator mf_ciup_begin, mf_ciup_end;
  op_iterator mf_cidn_begin, mf_cidn_end;
  op_iterator mf_niup_begin, mf_niup_end;
  op_iterator mf_nidn_begin, mf_nidn_end;
  op_iterator bsinglet_hop_begin, bsinglet_hop_end;

  // helper functions
  int set_parameter(const param& p, const std::string& name, const input::Parameters& inputs);
  //int set_mf_parameter(const mfparam& p, const std::string& name, const input::Parameters& inputs);
  MatrixElement product(const Complex& factor, const param& pname) const;
  MatrixElement product(const double& factor, const param& pname) const;
  //MatrixElement product(const double& factor, const param& pname);
  int add_operator(const op_id& quantum_op, const std::vector<MatrixElement>& MatrixElements, const unsigned& ntypes);
  std::pair<spin, bool> operator_type(const op_id& quantum_op) const;
  int read_mf_order(const input::Parameters& inputs);
  int read_pairing_symmetry(const input::Parameters& inputs);
  void finalize(void);
  void update_matrix_elements(void);
};


} // end namespace model

#endif