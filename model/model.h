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
#include <sstream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <stdexcept>
#include "../lattice/lattice.h"
#include "sitebasis.h"
#include "hamiltonian_term.h"

namespace model {

class Model : public std::vector<SiteTerm>,  public std::vector<BondTerm>
{
public:
  using siteterm_iterator = std::vector<SiteTerm>::const_iterator; 
  using bondterm_iterator = std::vector<BondTerm>::const_iterator; 
  Model() {}
  Model(const lattice::Lattice& lattice, const input::Parameters& inputs)
  { construct(lattice, inputs); }
  ~Model() {}
  int construct(const lattice::Lattice& lattice, const input::Parameters& inputs);

  unsigned add_sitebasis(SiteBasis& sitebasis, const unsigned& type=0);
  unsigned add_parameter(const std::string& pname, const double& defval, 
    const input::Parameters& inputs)
    { parms_[pname] = inputs.set_value(pname, defval); return parms_.size(); }
  void update_parameters(const input::Parameters& inputs);
  unsigned add_constant(const std::string& cname, const double& val) 
    { constants_.insert({cname, val}); return constants_.size(); }
  unsigned add_siteterm(const std::string& name, const CouplingConstant& cc,
    const std::string& op_expr, const std::string& site);
  unsigned add_bondterm(const std::string& name, const CouplingConstant& cc,
    const std::string& op_expr, const std::string& src, const std::string& tgt);

  const BasisDescriptor& basis(void) const { return basis_; }
  const SiteBasis& site_basis(const unsigned& site_type) const { return basis_.at(site_type); }
  unsigned sitebasis_dimension(const unsigned& site_type) const
  { return basis_.dimension(site_type); }
  const bool& has_siteterm(void) const { return has_siteterm_; }
  const bool& has_bondterm(void) const { return has_bondterm_; }
  const siteterm_iterator& siteterms_begin(void) const { return st_begin_; }
  const siteterm_iterator& siteterms_end(void) const { return st_end_; }
  const bondterm_iterator& bondterms_begin(void) const { return bt_begin_; }
  const bondterm_iterator& bondterms_end(void) const { return bt_end_; }
  std::pair<siteterm_iterator, siteterm_iterator> site_terms(void) const 
    { return std::make_pair(std::vector<SiteTerm>::cbegin(), std::vector<SiteTerm>::cend()); }
  std::pair<bondterm_iterator, bondterm_iterator> bond_terms(void) const 
    { return std::make_pair(std::vector<BondTerm>::cbegin(), std::vector<BondTerm>::cend()); }
  unsigned num_siteterms(void) const { return std::vector<SiteTerm>::size(); }
  unsigned num_bondterms(void) const { return std::vector<BondTerm>::size(); }
  unsigned num_total_terms(void) const 
    { return std::vector<SiteTerm>::size()+std::vector<BondTerm>::size(); }
  void get_term_names(std::vector<std::string>& term_names) const;
  std::ostream& print_info(std::ostream& os) const { return os << info_str_.str(); }

  //const SiteTerm& siteterm(const unsigned& i) const { return siteterms_[i]; };
  //const BondTerm& bondterm(const unsigned& i) const { return bondterms_[i]; };
private:
  int define_model(const lattice::Lattice& lattice, const input::Parameters& inputs);
  void finalize(const lattice::Lattice& lattice);
  void set_info_string(const lattice::Lattice& lattice); 

  std::string model_name;
  BasisDescriptor basis_;
  std::map<unsigned, unsigned> sitetypes_map_;
  std::map<unsigned, unsigned> bondtypes_map_;
  bool has_siteterm_{false};
  bool has_bondterm_{false};
  siteterm_iterator st_begin_;
  siteterm_iterator st_end_;
  bondterm_iterator bt_begin_;
  bondterm_iterator bt_end_;
  //std::vector<SiteTermDescriptor> siteterms_;
  //std::vector<BondTermDescriptor> bondterms_;
  ModelParams parms_;
  ModelParams constants_;

  std::ostringstream info_str_;
};


} // end namespace model

#endif