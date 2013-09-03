#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;

namespace G2G {
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
  cout << "iteration: " << t.total << endl;
  cout << "rmm: " << t.rmm << " density: " << t.density << " pot: " << t.pot << " forces: " << t.forces << " resto: " << t.resto << " functions: " << t.functions << endl;
  return io;
}

/********************
 * PointGroup
 ********************/

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          rmm_input(ii, jj) = (scalar_type)fortran_vars.rmm_input_ndens1.data[big_index];
          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          fortran_vars.rmm_output(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(), std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

template<class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
	double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
	double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
	delta = abs(x-x1);
	x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}

/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/

template class PointGroup<double>;
template class PointGroup<float>;
}