/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

std::pair<spin, bool> Model::operator_type(const op_id& op) const
{
  spin sigma;
  bool single_particle_op = true;
  switch (op) {
    case op_id::upspin_hop: sigma=spin::UP; single_particle_op=true; break;
    case op_id::dnspin_hop: sigma=spin::DOWN; single_particle_op=true; break;
    case op_id::ni_up: sigma=spin::UP; single_particle_op=true; break;
    case op_id::ni_dn: sigma=spin::DOWN; single_particle_op=true; break;
    case op_id::ni: sigma=spin::UD; single_particle_op=true; break;
    case op_id::Hubbard: sigma=spin::UD; single_particle_op=false; break;
    case op_id::Exchange: sigma=spin::UD; single_particle_op=false; break;
    case op_id::bondsinglet_hop: sigma=spin::UD; single_particle_op=false; break;
    default: throw std::range_error("*error: operator_type: unknown operator");
  }
  return std::pair<spin,bool>(sigma, single_particle_op);
}

int Model::set_parameter(const param& p, const std::string& name, const input::Parameters& inputs)
{
  parameters[p] = inputs.set_value(name, 0.0);
  return 0;
}

/*int Model::set_mf_parameter(const mfparam& p, const std::string& name, const input::Parameters& inputs)
{
  mf_parameters[p] = inputs.set_value(name, 0.0);
  return 0;
}*/

MatrixElement Model::product(const Complex& factor, const param& p) const
{
  return {Complex(0.0), factor, p};
}

MatrixElement Model::product(const double& factor, const param& p) const
{
  return {Complex(0.0), Complex(factor), p};
}

int Model::add_operator(const op_id& op, const std::vector<MatrixElement>& MatrixElements, const unsigned& ntypes)
{
  if (ntypes < 1) throw std::range_error("*error: add_operator: No of matrix element types must be > 0");
  std::vector<MatrixElement> matrix_elems(ntypes);
  for (unsigned i=0; i<ntypes; ++i) matrix_elems[i] = MatrixElements[i];
  spin sigma; bool single_particle_op;
  boost::tie(sigma, single_particle_op) = operator_type(op);
  operators.insert({op, op_value{spin::UP, matrix_elems, true}});
  return 0;
}

void Model::finalize(void) {
  update_matrix_elements();
  std::tie(mf_ciup_begin, mf_ciup_end) = mf_operators.equal_range(op_id::upspin_hop);
  std::tie(mf_cidn_begin, mf_cidn_end) = mf_operators.equal_range(op_id::dnspin_hop);
  std::tie(mf_niup_begin, mf_niup_end) = mf_operators.equal_range(op_id::ni_up);
  std::tie(mf_nidn_begin, mf_nidn_end) = mf_operators.equal_range(op_id::ni_dn);
  std::tie(bsinglet_hop_begin, bsinglet_hop_end) = mf_operators.equal_range(op_id::bondsinglet_hop);
}

void Model::update_matrix_elements(void)
{
  /*
  * For each operator,  matrix_element[type] = prefactor * model_parameter_value
  */
  for (auto it=operators.begin(); it!=operators.end(); ++it) {
    for (unsigned i = 0; i < it->second.matrix_elem.size(); ++i) {
      auto p = parameters.find(it->second.matrix_elem[i].pname);
      if (p != parameters.end()) {
        it->second.matrix_elem[i].value = it->second.matrix_elem[i].factor * p->second;
      }
      else {
        throw std::logic_error("*error: modellibrary: model parameter required by operator not set");
      }
    }
  }
  // update mean-field operators
  for (auto it=mf_operators.begin(); it!=mf_operators.end(); ++it) {
    for (unsigned i = 0; i < it->second.matrix_elem.size(); ++i) {
      auto p = parameters.find(it->second.matrix_elem[i].pname);
      if (p != parameters.end()) {
        it->second.matrix_elem[i].value = it->second.matrix_elem[i].factor * p->second;
      }
      else {
        throw std::logic_error("*error: modellibrary: model parameter required by operator not set");
      }
    }
  }
}



} // end namespace model
