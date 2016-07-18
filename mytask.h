/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 12:08:30
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-15 15:02:11
*----------------------------------------------------------------------------*/
#ifndef MYTASK_H
#define MYTASK_H

class MyTask : public scheduler::Task
{
public:
  MyTask() {};
  ~MyTask() {};
  void start(input::Parameters& parms); 
  void run() {}; 
  void dostep() {}; 
  void halt() {}; 
  //void run_test_code(const input::Parameters& parms); 
private:
  //int nsites;
};

#endif
