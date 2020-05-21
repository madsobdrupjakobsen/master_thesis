//
// Created by Niclas Laursen Brok on 2020-02-14.
//

#include <pybind11/pybind11.h>
#include "switching-times.hpp"

namespace SwitchingTimes {

    double cexp(double x, double cap) {
        double _expx = std::exp(x);
        if(x > cap) { _expx = std::exp(cap); };
        return _expx;
    };

    CppAD::AD<double> cexp(CppAD::AD<double> x, double cap) {
        CppAD::AD<double> _cap = CppAD::AD<double>(cap);
        return CppAD::CondExpGt(x, _cap, CppAD::exp(_cap), CppAD::exp(x));
    };

}