/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observable_operator.h"
#include "../expression/expression.h"

namespace mc {

SiteObsOperator::SiteObsOperator(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site)
{
  init(basis, op_expr, site);
}

SiteObsOperator::SiteObsOperator(const model::BasisDescriptor& basis, 
  const op_sitetypes& sitetypes, const std::string& op_expr, const std::string& site)
{
  init(basis, sitetypes, op_expr, site);
}

void SiteObsOperator::init(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site)
{
  // sitetype not specied, acts on all sites
  for (unsigned type=0; type<basis.size(); ++type) {
    push_back(model::SiteOperator(op_expr,site));
    back().build_matrix(basis.at(type));
  }
}

void SiteObsOperator::init(const model::BasisDescriptor& basis, 
  const op_sitetypes& sitetypes, const std::string& op_expr, const std::string& site)
{
  // construct the site operators
  for (unsigned type=0; type<basis.size(); ++type) {
    if (sitetypes.find(type) != sitetypes.end()) {
      push_back(model::SiteOperator(op_expr,site));
    }
    else push_back(model::SiteOperator());
    back().build_matrix(basis.at(type));
  }
}

} // end namespace mc

