//////////////////////////////////////////////////////////////////////////////////////////
//
//  Here are a handful of c++ functions which I created, simply to replace some fortran
//  functions which 
//
//
//////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <random>
#include <stdio.h>

#include "fortran_functs.hh"

using namespace std;


//__________________________________________________________________________________________
double fortran_functs::fcn_grnd() {
  
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(0., 1.);
  
  return dis(gen); 
}

//__________________________________________________________________________________________
double fortran_functs::fcn_gauss(double *std_dev) {
  
  random_device rd;
  mt19937 gen(rd());
  normal_distribution<double> dis(0., 1.);
  
  //  printf("std-dev passed: % 0.5e\n", *std_dev); 
  
  return dis(gen)*(*std_dev); 
}
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________

extern "C" {
  double grnd_(void) {
    return fortran_functs::fcn_grnd();
  }

  double gauss1_(double *std_dev) {
    return fortran_functs::fcn_gauss(std_dev);
  }
}
