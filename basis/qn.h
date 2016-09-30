/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef QN_H
#define QN_H

#include <iostream>
#include <string>
#include <vector>

namespace basis {

class QN_set : public std::vector<short> 
{
  private:
    static std::vector<std::string> qn_name;
    static int fermionic;
    static int qn_mask;
    static size_t QN_LAST;

    short qn0;
    short qn1;
    short qn2;
    short qn3;
public:
  QN_set();
  ~QN_set();
};

} // end namespace basis

#endif
