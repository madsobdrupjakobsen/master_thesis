//
// Created by Niclass Laursen Brok on 2020-02-15.
// Modified by Mads Obdrup Jakobsen on 2020-03-11.
//

#include "switching-times.hpp"

namespace SwitchingTimes {
    template<typename scalar>
    void Plant::model(const vector<scalar> &x, vector<scalar> &dxdt,
            const double t,
            const vector<scalar> &p_dynamic, const vector<scalar> &p_opt, const vector<double> &p_const) {
        /*
         * Extract dynamical parameters
         */
        Eigen::Map<const vector<scalar>> dap(p_dynamic.data(), 49);      // Day-ahead price
        Eigen::Map<const vector<scalar>> dat(p_dynamic.data() + 48, 49); // Day-ahead times
        /*
         * Fill model regime activation
         * p_opt = (ON-vec; OFF-vec)
         */
        size_t n_opt = p_opt.size() / 2;
        Eigen::Map<const vector<scalar>> on(p_opt.data(), n_opt);
        Eigen::Map<const vector<scalar>> off(p_opt.data() + n_opt, n_opt);
        scalar model_regime = 0.;
        for(int k = 0; k < n_opt; ++k) {
            model_regime += 1. / ((1. + cexp(-p_const(8) * (t - on(k)), 15.)) *
                                  (1. + cexp( p_const(8) * (t - off(k)), 15.)));
        };
        /*
         * Fill day-ahead price activation
         */
        scalar day_ahead_price = 0.;
        for(int k = 0; k < 48; ++k) {
            day_ahead_price += dap(k) / ((1. + cexp(-p_const(7) * (t - dat(k)), 15.)) *
                                         (1. + cexp( p_const(7) * (t - dat(k + 1)), 15.)));
        };
        /*
         * Compute dynamics
         */
        dxdt(0) =  model_regime * p_const(0) * (p_const(1) - x(0)) + 
                (1-model_regime) * p_const(2) * (p_const(3) - x(0));            // Process dynamic
        
        dxdt(1) = day_ahead_price * (p_const(4) * x(0) + 
                                    model_regime * p_const(5) + 
                                    (1-model_regime) * p_const(6)) * 1./60.;                         // Cost
    };
    template <typename scalar>
    scalar Plant::objective(const vector<scalar> &x,
            const vector<scalar> &p_dynamic, const vector<scalar> &p_opt, const vector<double> &p_const) {
        
        return (x(1));
    };

/*
 * Instantiate templates -> needed for proper linking at compile time
 */
template void Plant::model(const vector<double> &x, vector<double> &dxdt,
                           const double t,
                           const vector<double> &p_dynamic, const vector<double> &p_opt, const vector<double> &p_const);
template void Plant::model(const vector<ad_double> &x, vector<ad_double> &dxdt,
                           const double t,
                           const vector<ad_double> &p_dynamic, const vector<ad_double> &p_opt,
                           const vector<double> &p_const);
template double Plant::objective(const vector<double> &x,
                                 const vector<double> &p_dynamic, const vector<double> &p_opt,
                                 const vector<double> &p_const);
template ad_double Plant::objective(const vector<ad_double> &x,
                                    const vector<ad_double> &p_dynamic, const vector<ad_double> &p_opt,
                                    const vector<double> &p_const);

}