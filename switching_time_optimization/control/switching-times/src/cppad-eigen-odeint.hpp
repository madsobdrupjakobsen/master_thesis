//
// Created by Niclas Laursen Brok on 2020-02-14.
//

#ifndef SWITCHINGTIMES_CPPAD_EIGEN_ODEINT_HPP
#define SWITCHINGTIMES_CPPAD_EIGEN_ODEINT_HPP

#include <pybind11/eigen.h>
#include <cppad/cppad.hpp>
#include <Eigen/Dense>

// only for steppers with error control
namespace boost { namespace numeric { namespace odeint {
            template<int S1, int S2, int O, int M1, int M2>
            struct vector_space_norm_inf< Eigen::Matrix<CppAD::AD<double>, S1, S2, O, M1, M2> > {
            typedef double result_type;
            double operator()(const Eigen::Matrix<CppAD::AD<double>, S1, S2, O, M1, M2> &p) const {
                using std::max;
                using std::abs;
                result_type _out = abs(CppAD::Value(p(0)));
                for (int k = 1; k < p.size(); ++k) {
                    _out = max(_out, abs(CppAD::Value(p(k))));
                };
                return _out;
            }
        };
    } } }

#endif //SWITCHINGTIMES_CPPAD_EIGEN_ODEINT_HPP
