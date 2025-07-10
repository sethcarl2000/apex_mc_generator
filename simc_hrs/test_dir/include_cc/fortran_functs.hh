#ifndef FORTRAN_FUNCTS_HH
#define FORTRAN_FUNCTS_HH

#include <iostream>
#include <random>
#include <stdio.h>

//generates a random number on the interval [0,1]
namespace fortran_functs {

  double fcn_grnd();

  double fcn_gauss(double *std_dev); 
};

#endif
