/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "model.h"

namespace model {

unsigned Model::add_sitebasis(SiteBasis& sitebasis, const unsigned& type)
{
  auto it=sitetypes_map_.find(type);
  if (it==sitetypes_map_.end()) 
    throw std::range_error("Model::add_sitebasis: specified 'site type' not found");
  unsigned mapped_type = it->second;
  return basis_.add_sitebasis(mapped_type,sitebasis); 
}

unsigned Model::add_siteterm(const std::string& name, const CouplingConstant& cc,
  const std::string& op_expr, const std::string& site)
{
  // remap site type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  for (auto it=cc.begin(); it!=cc.end(); ++it) {
    unsigned sitetype = it->first;
    auto it2=sitetypes_map_.find(sitetype);
    if (it2!=sitetypes_map_.end()) {
      unsigned mapped_type = it2->second;
      cc_remapped.insert({mapped_type, it->second});
    }
    else throw std::range_error("Model::add_siteterm: non-existent 'site type' specified");
  }
  unsigned num_sitetypes = sitetypes_map_.size();
  this->std::vector<SiteTerm>::push_back(SiteTerm(name, cc_remapped, op_expr, site, num_sitetypes));
  return this->std::vector<SiteTerm>::size();
}

unsigned Model::add_bondterm(const std::string& name, const CouplingConstant& cc,
  const std::string& op_expr, const std::string& src, const std::string& tgt)
{
  // remap bond type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  for (auto it=cc.begin(); it!=cc.end(); ++it) {
    unsigned bondtype = it->first;
    auto it2=bondtypes_map_.find(bondtype);
    if (it2!=bondtypes_map_.end()) {
      unsigned mapped_type = it2->second;
      cc_remapped.insert({mapped_type, it->second});
    }
    else throw std::range_error("Model::add_bondterm: non-existent 'site type' specified");
  }
  unsigned num_bondtypes = bondtypes_map_.size();
  std::vector<BondTerm>::push_back(BondTerm(name, cc_remapped, op_expr, src, tgt, num_bondtypes));
  return std::vector<BondTerm>::size();
}

void Model::finalize(const lattice::Lattice& L)
{
  // finalize the site terms
  for (auto it=std::vector<SiteTerm>::begin(); it!=std::vector<SiteTerm>::end(); ++it) {
    it->build_matrix(basis_); 
    it->eval_coupling_constant(constants_, parms_); 
  }
  has_siteterm_ = (std::vector<SiteTerm>::size()>0);
  st_begin_ = std::vector<SiteTerm>::cbegin();
  st_end_ = std::vector<SiteTerm>::cend();
  // finalize the bond terms
  // map each 'bondtype' to its 'src' & 'tgt' types
  BondTerm::BondSiteMap bond_site_map;
  for (unsigned i=0; i<L.num_unitcell_bonds(); ++i) {
    lattice::Bond b = L.unitcell_bond(i);
    lattice::Site src = L.unitcell_site(b.src_id());
    lattice::Site tgt = L.unitcell_site(b.tgt_id());
    bond_site_map.insert({b.type(), std::make_pair(src.type(), tgt.type())});
  }
  for (auto it=std::vector<BondTerm>::begin(); it!=std::vector<BondTerm>::end(); ++it) {
    it->build_matrix(basis_, bond_site_map); 
    it->eval_coupling_constant(constants_, parms_); 
  }
  has_bondterm_ = (std::vector<BondTerm>::size()>0);
  bt_begin_ = std::vector<BondTerm>::cbegin();
  bt_end_ = std::vector<BondTerm>::cend();

  set_info_string(L);
}

void Model::update_parameters(const input::Parameters& inputs)
{
  // update the parameter values
  for (auto& p : parms_) p.second = inputs.set_value(p.first, p.second);
  // update the model term couping constants
  for (auto it=std::vector<SiteTerm>::begin(); it!=std::vector<SiteTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  for (auto it=std::vector<BondTerm>::begin(); it!=std::vector<BondTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
}

void Model::get_term_names(std::vector<std::string>& term_names) const
{
  term_names.clear();
  for (auto it=std::vector<BondTerm>::cbegin(); it!= std::vector<BondTerm>::cend(); ++it) 
    term_names.push_back(it->name());
  for (auto it=std::vector<SiteTerm>::cbegin(); it!= std::vector<SiteTerm>::cend(); ++it) 
    term_names.push_back(it->name());
}

void Model::set_info_string(const lattice::Lattice& L) 
{
  info_str_.clear();
  info_str_ << "# Lattice: " << L.name() << " (";
  info_str_ << "Size="<<L.size1()<<"x"<<L.size2()<<"x"<< L.size3()<<", ";
  info_str_ << "Sites/unitcell="<<L.num_basis_sites()<<", ";
  info_str_ << "Boundary="<<static_cast<int>(L.bc1_periodicity()) << "-"; 
  info_str_ << static_cast<int>(L.bc2_periodicity()) << "-";
  info_str_ << static_cast<int>(L.bc3_periodicity()) << ")\n";
  info_str_ << "# No of sites = " << L.num_sites() << "\n";
  info_str_ << "# Model: " << model_name << "\n";
  info_str_.precision(6);
  info_str_.setf(std::ios_base::fixed);
  for (const auto& p : parms_) 
    info_str_ << "# " << p.first << " = " << p.second << "\n";
}





} // end namespace model
