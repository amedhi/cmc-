/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "hamiltonian_term.h"

namespace model {

/*--------------------------CouplingConstant--------------------*/
const int CouplingConstant::global_type = -1;

CouplingConstant::CouplingConstant(const std::string& expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1;
  valid_ = true;
}

std::pair<CouplingConstant::iterator, bool> CouplingConstant::insert(const value_type& val) 
{
  std::pair<iterator, bool> res = super_type::insert(val);
  num_types_ = super_type::size();
  valid_ = true;
  return res;
}

CouplingConstant& CouplingConstant::operator=(const std::string expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

CouplingConstant::CouplingConstant(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  create(type0, type1, type2, type3, type4, type5); 
}

void CouplingConstant::clear(void)
{
  super_type::clear();
  num_types_ = 0;
  valid_ = false;
}

void CouplingConstant::create(const unsigned& num_types) 
{
  super_type::clear();
  num_types_ = num_types;
  valid_ = false;
}

void CouplingConstant::create(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  super_type::insert(type0);
  num_types_=1;
  if (type1.second != "_null_") {
    super_type::insert(type1); num_types_++;
  }
  if (type2.second != "_null_") {
    super_type::insert(type2); num_types_++;
  }
  if (type3.second != "_null_") {
    super_type::insert(type3); num_types_++;
  }
  if (type4.second != "_null_") {
    super_type::insert(type4); num_types_++;
  }
  if (type5.second != "_null_") {
    super_type::insert(type5); num_types_++;
  }
  valid_=true;
}

void CouplingConstant::add_type(const unsigned& type, const std::string& expr)
{
  super_type::insert({type, expr});
  valid_ = (num_types_==static_cast<int>(size()));
}

void CouplingConstant::add_type(const value_type& val)
{
  super_type::insert(val);
  valid_ = (num_types_==static_cast<int>(size()));
}


//-----------------------SiteOperator---------------------------
SiteOperator::SiteOperator(const std::string& op_expr, const std::string& site) 
  : op_expr_(op_expr), site_(site)
{
  // strip the '(site)' substrings
  /*std::string s = "("+site_+")";
  while(true) {
    std::string::size_type pos = op_expr_.find(s);
    if (pos == std::string::npos) break;
    op_expr_.erase(pos,s.size());
  }*/
}

int SiteOperator::build_matrix(const SiteBasis& site_basis)
{
  matrix_.resize(site_basis.dimension());
  //std::cout << "expression = " << op_expr_ << std::endl;
  if (op_expr_.size()==0) {
    // this operator does not operate on this site
    for (unsigned i=0; i<site_basis.dimension(); ++i) matrix_(i) = 0.0;
    return 0;
  }

  // oprator expression: strip '(' & ')' characters around 'site'
  std::string expr_str(op_expr_);
  std::string s = "("+site_+")";
  std::string::size_type pos;
  while ((pos=expr_str.find(s)) != std::string::npos) expr_str.replace(pos,s.length(),site_);
  //while ((pos=expr_str.find('(')) != std::string::npos) expr_str.erase(pos,1);
  //while ((pos=expr_str.find(')')) != std::string::npos) expr_str.erase(pos,1);
  //std::cout << "expression = " << expr_str << std::endl;

  // expression evaluator
  expr::Expression op_expr; 
  expr::Expression::variables vars;

  unsigned num_operator = site_basis.num_operator();
  for (unsigned i=0; i<site_basis.dimension(); ++i) {
    for (unsigned n=0; n<num_operator; ++n) {
      QuantumOperator op = site_basis.quantum_operator(n);
      std::string vname = op.name()+site_;
      vars[vname] = static_cast<double>(site_basis.apply(op.id(), i));
      //expr.set_variable(op.name(), m);
    }
    matrix_(i) = op_expr.evaluate(expr_str, vars);
    //std::cout << "state, val =" << i << ", " << static_cast<int>(res) << std::endl;
  }
  return 0;
}

// when default variable values are set
int SiteOperator::build_matrix(const SiteBasis& site_basis, expr::Expression::variables& vars)
{
  matrix_.resize(site_basis.dimension());
  //std::cout << "expression = " << op_expr_ << std::endl;
  if (op_expr_.size()==0) {
    // this operator does not operate on this site
    for (unsigned i=0; i<site_basis.dimension(); ++i) matrix_(i) = 0.0;
    return 0;
  }

  // oprator expression: strip '(' & ')' characters around 'site'
  std::string expr_str(op_expr_);
  std::string s = "("+site_+")";
  std::string::size_type pos;
  while ((pos=expr_str.find(s)) != std::string::npos) expr_str.replace(pos,s.length(),site_);
  //while ((pos=expr_str.find('(')) != std::string::npos) expr_str.erase(pos,1);
  //while ((pos=expr_str.find(')')) != std::string::npos) expr_str.erase(pos,1);
  //std::cout << "expression = " << expr_str << std::endl;

  // expression evaluator
  expr::Expression op_expr; 

  unsigned num_operator = site_basis.num_operator();
  for (unsigned i=0; i<site_basis.dimension(); ++i) {
    for (unsigned n=0; n<num_operator; ++n) {
      QuantumOperator op = site_basis.quantum_operator(n);
      std::string vname = op.name()+site_;
      vars[vname] = static_cast<double>(site_basis.apply(op.id(), i));
      //expr.set_variable(op.name(), m);
    }
    matrix_(i) = op_expr.evaluate(expr_str, vars);
    //std::cout << "state, val =" << i << ", " << static_cast<int>(res) << std::endl;
  }
  return 0;
}

//-----------------------SiteOperatorTerm-------------------------
/*SiteOperatorTerm::SiteOperatorTerm(const std::string& name, const double& coupling, 
  const std::string& site, const std::string& op_expr)
  : SiteOperator(site, op_expr), name_{name}, coupling_{coupling}
{
}
*/

SiteOperatorTerm::SiteOperatorTerm(const std::string& name, const std::string& cc_expr, 
  const std::string& op_expr, const std::string& site)
  : SiteOperator(op_expr, site), name_{name}, cc_expr_{cc_expr}
{
}

int SiteOperatorTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  if (cc_expr_.size()==0) {
    cc_value_ = 0.0;
    return 0;
  }

  // expression evaluator
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& c : cvals) vars[c.first] = c.second;
  for (const auto& p : pvals) vars[p.first] = p.second;
  try { 
    cc_value_ = expr.evaluate(cc_expr_, vars); 
  }
  catch (std::exception& e) 
  { 
    std::string msg = "SiteOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
  return 0;
}

//-----------------------SiteTerm-------------------------
/*SiteTerm::SiteTerm(const std::string& name, const TermParameter& coupling, 
    const std::string& site, const std::string& op_expr)
{
  if (!coupling.valid()) throw std::invalid_argument("SiteTerm::Invalid TermParameter");
  for (const auto& p : coupling) {
    std::string term_name = name + std::to_string(p.first);
    insert({p.first, SiteOperatorTerm(term_name,p.second,site,op_expr)});
  }
  name_ = name;
}*/

SiteTerm::SiteTerm(const std::string& name, const CouplingConstant& cc, 
    const std::string& op_expr, const std::string& site, const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("SiteTerm::Invalid CouplingConstant");
  name_ = name;
  size_ = size;

  // if the 'cc' is implicitly defined for all site types 
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    std::string cc_expr = cc.begin()->second;
    for (unsigned i=0; i<size_; ++i) {
      std::string term_name = name + std::to_string(i);
      this->operator[](i) = SiteOperatorTerm(term_name,cc_expr,op_expr,site);
    }
  } 
  else {
    // set "SiteOperatorTerm"-s only for those site types for which 'cc' is 
    // explicitly defined
    for (const auto& p : cc) {
      std::string term_name = name + std::to_string(p.first);
      std::string cc_expr = p.second;
      //std::cout << p.first << " " << p.second << "\n";
      this->operator[](p.first) = SiteOperatorTerm(term_name,cc_expr,op_expr,site);
      //insert({p.first, SiteOperatorTerm(term_name,p.second,op_expr,site)});
    }
  }
}

void SiteTerm::build_matrix(const BasisDescriptor& basis) 
{
  for (unsigned i=0; i<size_; ++i)
    this->operator[](i).build_matrix(basis.site_basis(i));
  //for (auto& elem : *this) {
  //  unsigned type = elem.first;
  //  elem.second.build_matrix(basis.site_basis(type));
  //}
}

void SiteTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  for (auto& term : *this) term.eval_coupling_constant(cvals, pvals);
  //for (auto& elem : *this) {
  //  elem.second.eval_coupling_constant(pvals);
  //}
}

const double& SiteTerm::matrix_element(const unsigned& site_type, const unsigned& state_idx) const
{
  return this->operator[](site_type).matrix_element(state_idx);
  //const_iterator it=find(site_type);
  //if (it != end()) return it->second.matrix_element(state_idx);
  //else return null_element_;
  //return at(site_type).matrix_element(state_idx);
} 

const double& SiteTerm::coupling(const unsigned& site_type) const
{
  return this->operator[](site_type).coupling();
  //const_iterator it=find(site_type);
  //if (it != end()) return it->second.coupling();
  //else return null_coupling_;
}


//-----------------------BondOperator---------------------------
BondOperator::BondOperator(const std::string& op_expr, const std::string& source, 
  const std::string& target) : op_expr_(op_expr), source_(source), target_(target)
{
}

int BondOperator::build_matrix(const SiteBasis& source_basis, const SiteBasis& target_basis)
{

  unsigned M = source_basis.dimension();
  unsigned N = target_basis.dimension();
  //std::cout <<"M= "<<M<<" N= "<<N<<"\n"; 
  matrix_.resize(M,N);
  if (op_expr_.size()==0) {
    // this operator does not operate on this bond
    for (unsigned i=0; i<M; ++i) 
      for (unsigned j=0; j<N; ++j) 
        matrix_(i,j) = 0.0;
    return 0;
  }

  // oprator expression: strip '(' & ')' characters
  std::string expr_str(op_expr_);
  std::string::size_type pos;
  std::string s = "("+source_+")";
  while ((pos=expr_str.find(s)) != std::string::npos) expr_str.replace(pos,s.length(),source_);
  s = "("+target_+")";
  while ((pos=expr_str.find(s)) != std::string::npos) expr_str.replace(pos,s.length(),target_);
  //while ((pos=expr_str.find('(')) != std::string::npos) expr_str.erase(pos,1);
  //while ((pos=expr_str.find(')')) != std::string::npos) expr_str.erase(pos,1);
  //std::cout << "expression = " << expr_str << std::endl;

  // expression evaluator
  expr::Expression op_expr; //(expr_str);
  expr::Expression::variables vars;

  for (unsigned i=0; i<M; ++i) {
    for (unsigned m=0; m<source_basis.num_operator(); ++m) {
      QuantumOperator op = source_basis.quantum_operator(m);
      std::string vname = op.name()+source_;
      vars[vname] = static_cast<double>(source_basis.apply(op.id(), i));
    }
    for (unsigned j=0; j<N; ++j) {
      for (unsigned n=0; n<target_basis.num_operator(); ++n) {
        QuantumOperator op = target_basis.quantum_operator(n);
        std::string vname = op.name()+target_;
        vars[vname] = static_cast<double>(target_basis.apply(op.id(), j));
      }
      matrix_(i,j) = op_expr.evaluate(expr_str, vars);
    }
  }
  return 0;
}

//-----------------------BondOperatorTerm-------------------------
BondOperatorTerm::BondOperatorTerm(const std::string& name, const std::string& cc_expr, 
  const std::string& op_expr, const std::string& source, const std::string& target)
  : BondOperator(op_expr, source, target), name_{name}, cc_expr_{cc_expr}, cc_value_{0.0}
{
}

int BondOperatorTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  if (cc_expr_.size()==0) {
    cc_value_ = 0.0;
    return 0;
  }
  // expression evaluator
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& c : cvals) vars[c.first] = c.second;
  for (const auto& p : pvals) vars[p.first] = p.second;
  try { 
    cc_value_ = expr.evaluate(cc_expr_, vars); 
    //std::cout << "bterm: " << name_ << " = " << cc_value_ << "\n";
  }
  catch (std::exception& e) 
  { 
    std::string msg = "BondOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
  return 0;
}


//-----------------------BondTerm-------------------------
BondTerm::BondTerm(const std::string& name, const CouplingConstant& cc, 
    const std::string& op_expr, const std::string& src, const std::string& tgt, 
    const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("BondTerm:: Invalid CouplingConstant");
  name_ = name;
  size_ = size;

  // if the 'cc' is implicitly defined for all bond types 
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    std::string cc_expr = cc.begin()->second;
    for (unsigned i=0; i<size_; ++i) {
      std::string term_name = name + std::to_string(i);
      this->operator[](i) = BondOperatorTerm(term_name,cc_expr,op_expr,src,tgt);
    }
  } 
  else {
    // set "BondOperatorTerm"-s only for those bond types for which 'cc' is 
    // defined explicitly
    for (const auto& p : cc) {
      std::string term_name = name + std::to_string(p.first);
      std::string cc_expr = p.second;
      this->operator[](p.first) = BondOperatorTerm(term_name,cc_expr,op_expr,src,tgt);
      //insert({p.first, BondOperatorTerm(term_name,p.second,op_expr,src,tgt)});
    }
  }
}

void BondTerm::build_matrix(const BasisDescriptor& basis, const BondSiteMap& bondtypes)
{
  // build matrix for all the bond types in the lattice 
  // (for types with undefined 'cc' matrices are zero)
  unsigned src_type, tgt_type;
  for (unsigned i=0; i<size_; ++i) {
    std::tie(src_type, tgt_type) = bondtypes.at(i);
    //std::cout << size_ <<" "<<i<<" "<<src_type<<" "<<tgt_type<<"\n"; 
    this->operator[](i).build_matrix(basis.site_basis(src_type), basis.site_basis(tgt_type));
  }
  /*for (auto& elem : *this) {
    unsigned btype = elem.first;
    BondSiteMap::const_iterator it = bondtypes.find(btype);
    if (it != bondtypes.cend()) {
      std::tie(src_type, tgt_type) = it->second;
      elem.second.build_matrix(basis.site_basis(src_type), basis.site_basis(tgt_type));
    }
    else
      throw std::range_error("BondTerm::build_matrix: undefined bond type");
  }*/
}

void BondTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  for (auto& term : *this) term.eval_coupling_constant(cvals, pvals);
  //std::cout << "------hi--------\n"; //abort();
  /*for (auto& elem : *this) {
    elem.second.eval_coupling_constant(pvals);
  }*/
}

const double& BondTerm::matrix_element(const unsigned& bond_type, const unsigned& src_idx,
  const unsigned& tgt_idx) const
{
  return this->operator[](bond_type).matrix_element(src_idx, tgt_idx);
  /*const_iterator it=find(bond_type);
  if (it != end()) return it->second.matrix_element(src_idx, tgt_idx);
  else return null_element_;
  */
} 

const double& BondTerm::coupling(const unsigned& bond_type) const
{
  return this->operator[](bond_type).coupling();
  /*const_iterator it=find(bond_type);
  if (it != end()) return it->second.coupling();
  else return null_coupling_;
  */
}



} // end namespace model
