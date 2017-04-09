/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-09 01:45:56
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Model::define_model(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  //int info;
  //unsigned ntypes;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  unsigned id;
  unsigned sitetype, change, type, src_type, tgt_type;
  std::string name, matrixelem, op, qn, site, src, tgt, fact;
  SiteBasis site_basis;
  BasisDescriptor basis;
  QuantumNumber::value_type min, max, step;
  CouplingConstant cc;
  //QuantumOperator q_op0, q_op1, q_op2;

  // define the models 
  model_name = inputs.set_value("model", "ISING");
  boost::to_upper(model_name);

  if (model_name == "ISING") {
    mid = model_id::ISING;
    // site basis
    site_basis.clear();
    site_basis.add_qn(qn="S", min=-1, max=1, step=2);
    site_basis.add_operator(op="S", matrixelem="S", qn="S", change=0);
    add_sitebasis(site_basis);
    // model parameters
    add_parameter(name="J", defval=1.0, inputs);
    // bond operator terms
    add_bondterm(name="Exchange", cc="-J", op="S(i)*S(j)", src="i", tgt="j");
  }

  else if (model_name == "POTTS") {
    mid = model_id::POTTS;

    // q-value of the model
    int q = inputs.set_value("q", 2);
    if (q < 2) throw std::range_error("modellibrary: POTTS model 'q' must be >= 2");
    add_constant("q", q);

    // site basis
    site_basis.clear();
    site_basis.add_qn(qn="S", min=0, max=(q-1), step=1);
    site_basis.add_operator(op="S", matrixelem="S", qn="S");
    add_sitebasis(site_basis);
    // model parameters
    add_parameter(name="J", defval=1.0, inputs);
    add_parameter(name="H", defval=0.0, inputs);
    // site operator term
    // magnetic field along +z direction 
    add_siteterm("H_field", cc="-H", op="cron(S(i),0)", site="i");
    // bond operator terms
    add_bondterm(name="Exchange", cc="-J", op="cron(S(i),S(j))", src="i", tgt="j");
  }

  else if (model_name == "BEG") {
    mid = model_id::BEG;
    //site_basis
    site_basis.clear();
    site_basis.add_qn(qn="sigma", min=-1, max=1, step=1);
    site_basis.add_operator(op="sigma", matrixelem="sigma", qn="sigma");
    add_sitebasis(site_basis);
    // model parameters
    add_parameter(name="K", defval=1.0, inputs);
    // bond operator term
    add_bondterm("Tetra", cc="-1.0", op="sigma(i)*sigma(j)", src="i", tgt="j");
    add_bondterm("Cubic", cc="-K", op="(1-sigma(i)*sigma(i))*(1-sigma(j)*sigma(j))", src="i", tgt="j");
    // site operator term
    // add_siteterm("Site", cc="4.0*K", op="sigma(i)*sigma(i)", src="i");
  }

  else if (model_name == "BEG_POTTS") {
    mid = model_id::BEG_POTTS;
    switch (lattice.id()) {
      /*------------- 'SIMPLE CUBIC' lattice--------------*/
      case lattice::lattice_id::SIMPLECUBIC:
        //site_basis
        site_basis.clear();
        site_basis.add_qn(qn="S", min=-2, max=2, step=1);
        site_basis.add_qn(qn="sigma", min=-1, max=1, step=1);
        site_basis.add_operator(op="S", matrixelem="S", qn="S");
        site_basis.add_operator(op="sigma", matrixelem="sigma", qn="sigma");
        add_sitebasis(sitetype=0, site_basis);

        // model parameters
        add_parameter(name="T", defval=1.0, inputs);
        add_parameter(name="kB", defval=1.0, inputs);
        add_parameter(name="J", defval=1.0, inputs);
        add_parameter(name="J_fm", defval=1.0, inputs);
        add_parameter(name="J_afm", defval=1.0, inputs);
        add_parameter(name="H", defval=0.0, inputs);
        add_parameter(name="U", defval=1.0, inputs);
        add_parameter(name="K", defval=1.0, inputs);
        add_parameter(name="T_afm", defval=1.0, inputs);
        // constants
        add_constant("g", 2.0);
        add_constant("muB", 0.057884);  // meV/Tesla
        add_constant("ln_p", std::log(2.0));  // ln(p)

        // this model has terms for impurity 'bond' 
        add_impurity_bond(id=0, type=1, src_type=0, tgt_type=0);
        //def_impurity_bondtype(type=1, src_type=0, tgt_type=0);

        // site operator term
        // magnetic field along +z direction (S=+2) 
        /*
        add_siteterm("H_field", cc="-g*muB*H/J", op="cron(S(i),2)", site="i");
        add_siteterm("sigma", cc="-ln_p*kB*T/J", op="1-sigma(i)*sigma(i)", site="i");

        // bond operator term
        cc = CouplingConstant({0, "-J_fm/J"}, {1, "-J_afm/J*min(1.0,(T/T_afm-1.0))"});
        add_bondterm("Potts", cc, op="cron(S(i),S(j))", src="i", tgt="j");
        add_bondterm("Tetra", cc="-1.0", op="sigma(i)*sigma(j)", src="i", tgt="j");
        add_bondterm("Cubic", cc="-K/J", op="(1-sigma(i)*sigma(i))*(1-sigma(j)*sigma(j))", src="i", tgt="j");
        add_bondterm("Interaction", cc="0.5*U/J", 
          op="cron(S(i),S(j))*((1-2*sigma(i)*sigma(i))*(1-2*sigma(j)*sigma(j))-1.0)", src="i", tgt="j");
        */
        add_siteterm("H_field", cc="-g*muB*H", op="cron(S(i),2)", site="i");
        add_siteterm("sigma", cc="-kB*T*ln_p", op="1.0-sigma(i)*sigma(i)", site="i");

        // bond operator term
        cc = CouplingConstant({0, "-J_fm"}, {1, "-J_afm*min(1.0,(T/T_afm-1.0))"});
        add_bondterm("Potts", cc, op="cron(S(i),S(j))", src="i", tgt="j");
        add_bondterm("Tetra", cc="-J", op="sigma(i)*sigma(j)", src="i", tgt="j");
        add_bondterm("Cubic", cc="-K", op="(1.0-sigma(i)*sigma(i))*(1.0-sigma(j)*sigma(j))", src="i", tgt="j");
        add_bondterm("Interaction", cc="0.5*U", 
          op="cron(S(i),S(j))*((1-2*sigma(i)*sigma(i))*(1-2*sigma(j)*sigma(j))-1.0)", src="i", tgt="j");
        break;

      /*------------- 'NIMNX' lattice--------------*/
      case lattice::lattice_id::SYS_NIMNX:
        //site_basis
        // X-atom sites: only structural degrees of freedom
        site_basis.clear();
        //site_basis.add_qn(qn="S", min=0, max=0, step=0);
        site_basis.add_qn(qn="sigma", min=-1, max=1, step=1);
        //site_basis.add_operator(op="S", matrixelem="S", qn="S");
        site_basis.add_operator(op="sigma", matrixelem="sigma", qn="sigma");
        add_sitebasis(sitetype=0, site_basis);
        // Mn-atom sites: structural & 5 state spin variables
        site_basis.clear();
        site_basis.add_qn(qn="S", min=-2, max=2, step=1);
        site_basis.add_qn(qn="sigma", min=-1, max=1, step=1);
        site_basis.add_operator(op="S", matrixelem="S", qn="S");
        site_basis.add_operator(op="sigma", matrixelem="sigma", qn="sigma");
        add_sitebasis(sitetype=1, site_basis);
        // Ni-atom sites: structural & 3 state spin variables
        site_basis.clear();
        site_basis.add_qn(qn="S", min=-1, max=1, step=1);
        site_basis.add_qn(qn="sigma", min=-1, max=1, step=1);
        site_basis.add_operator(op="S", matrixelem="S", qn="S");
        site_basis.add_operator(op="sigma", matrixelem="sigma", qn="sigma");
        add_sitebasis(sitetype=2, site_basis);

        // declare the impurity bond types 
        // replace 'type-1 & type-2' bonds
        add_impurity_bond(id=0, type=5, src_type=2, tgt_type=2);
        // replace 'type-3' bonds (with impurity site as source)
        add_impurity_bond(id=1, type=6, src_type=2, tgt_type=1);
        // replace 'type-3' bonds (with impurity site as target)
        add_impurity_bond(id=2, type=7, src_type=1, tgt_type=2);
        // replace 'type-4' bonds (with impurity site as source)
        add_impurity_bond(id=3, type=8, src_type=2, tgt_type=1);
        // replace 'type-4' bonds (with impurity site as target)
        add_impurity_bond(id=4, type=9, src_type=1, tgt_type=2);

        // model parameters
        add_parameter(name="T", defval=1.0, inputs);
        add_parameter(name="kB", defval=1.0, inputs);
        add_parameter(name="J", defval=1.0, inputs);
        add_parameter(name="K", defval=1.0, inputs);
        add_parameter(name="Jm_MnNi", defval=1.0, inputs);
        add_parameter(name="Jm_MnMn3", defval=1.0, inputs);
        add_parameter(name="Jm_MnMn4", defval=1.0, inputs);
        add_parameter(name="Jm_MnNiB", defval=1.0, inputs);
        add_parameter(name="U_MnNi", defval=1.0, inputs);
        add_parameter(name="U_MnMn3", defval=1.0, inputs);
        add_parameter(name="U_MnMn4", defval=1.0, inputs);
        add_parameter(name="U_MnNiB", defval=1.0, inputs);
        add_parameter(name="H", defval=0.0, inputs);
        // constants
        add_constant("g", 2.0);
        add_constant("muB", 0.057884);  // meV/Tesla
        add_constant("ln_p", std::log(2.0));  // ln(p)

        // site operator term
        cc = CouplingConstant({1, "-g*muB*H"}, {2, "-g*muB*H"});
        add_siteterm("H_field", cc, op="cron(S(i),0)", site="i");
        //add_siteterm("sigma", cc="-kB*T*ln_p", op="1.0-sigma(i)*sigma(i)", site="i");

        // bond operator term
        // magnetic exchange
        //cc = CouplingConstant({1,"-Jm_MnNi"}, {2,"-Jm_MnNi"}, {3,"-Jm_MnMn"},{4,"-Jm_MnMn4"});
        cc.create(6);
        cc.add_type(1, "-Jm_MnNi");
        cc.add_type(2, "-Jm_MnNi");
        cc.add_type(3, "-Jm_MnMn3");
        cc.add_type(4, "-Jm_MnMn4");
        cc.add_type(6, "-Jm_MnNiB");
        cc.add_type(7, "-Jm_MnNiB");
        add_bondterm("Potts", cc, op="cron(S(i),S(j))", src="i", tgt="j");

        // elastic J
        cc = CouplingConstant({0,"-J"}, {1,"-J"}, {2,"-J"}, {5,"-J"});
        add_bondterm("Tetra", cc="-J", op="sigma(i)*sigma(j)", src="i", tgt="j");
        // elastic K
        cc = CouplingConstant({0,"-K"}, {1,"-K"}, {2,"-K"}, {5,"-K"});
        add_bondterm("Cubic", cc="-K", op="(1.0-sigma(i)*sigma(i))*(1.0-sigma(j)*sigma(j))", src="i", tgt="j");
        // interaction
        //cc = CouplingConstant({1,"-U_MnNi"}, {2,"-U_MnNi"}, {3,"-U_MnMn"}, {4,"-U_MnMn4"});
        cc.create(6);
        cc.add_type(1, "-U_MnNi");
        cc.add_type(2, "-U_MnNi");
        cc.add_type(3, "-U_MnMn3");
        cc.add_type(4, "-U_MnMn4");
        cc.add_type(6, "-U_MnNiB");
        cc.add_type(7, "-U_MnNiB");
        add_bondterm("Interaction", cc, 
          op="cron(S(i),S(j))*sigma(i)*sigma(j)", src="i", tgt="j");
        //add_bondterm("MagElastic", cc="-K1", op="(1.0-sigma(i)*sigma(i))*(1.0-sigma(j)*sigma(j))", src="i", tgt="j");
        //add_bondterm("Interaction", cc="-U", 
        //  op="cron(S(i),S(j))*sigma(i)*sigma(i)*sigma(j)*sigma(j)", src="i", tgt="j");
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

int Model::construct(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  // reset
  parms_.clear();
  //operators.clear();
  // maps of site & bond type values (to contigous type values)
  sitetypes_map_ = lattice.sitetypes_map();
  bondtypes_map_ = lattice.bondtypes_map();
  // maps of a given bond type to the types of its target
  bond_sites_map_.clear();
  for (unsigned i=0; i<lattice.num_unitcell_bonds(); ++i) {
    lattice::Bond b = lattice.unitcell_bond(i);
    lattice::Site src = lattice.unitcell_site(b.src_id());
    lattice::Site tgt = lattice.unitcell_site(b.tgt_id());
    bond_sites_map_.insert({b.type(), std::make_pair(src.type(), tgt.type())});
    //std::cout << "bond_site_map = "<<b.type()<<" "<<src.type()<<" "<<tgt.type()<<"\n";
  }
  // impurity 
  impurity_bond_types_.clear();

  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model
