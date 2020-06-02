//
// Created by MadNiclass Laursen Brok on 2020-02-15.
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
            model_regime += 1. / ((1. + cexp(-p_const(12) * (t - on(k)), 15.)) *
                                  (1. + cexp( p_const(12) * (t - off(k)), 15.)));
        };
        /*
         * Fill day-ahead price activation
         */
        scalar day_ahead_price = 0.;
        for(int k = 0; k < 48; ++k) {
            day_ahead_price += dap(k) / ((1. + cexp(-p_const(11) * (t - dat(k)), 15.)) *
                                         (1. + cexp( p_const(11) * (t - dat(k + 1)), 15.)));
        };
        /*
         * Compute dynamics
         
        dxdt(0) =  x(1) / 889.190900140 ;

        dxdt(1) = model_regime * (-2 * p_const(5) * p_const(3) * (x(1)) - p_const(3)*p_const(3) * (889.190900140*x(0)-698.861205845*p_const(1))) +
                (1-model_regime)*(-2 * p_const(4) * p_const(2) * (x(1)) - p_const(2)*p_const(2) * (889.190900140*x(0)-961.996347735*p_const(0)));          // Process dynamic
        
        dxdt(2) =  day_ahead_price * (p_const(8) * 1./(1. + cexp(-(p_const(6) * (889.190900140*x(0) - 592.010123492*p_const(7))),15)) +  // Divide by 1000 to get into [0,1]
                                    model_regime * p_const(9) + 
                                    (1-model_regime) * p_const(10)) * 1./60.;                   // Cost
        */

        // dxdt[0] = x[1] / 889.190900140 
        // dxdt[1] = -2 * self.xi_MELT * self.omega_MELT * (x[1]) - self.omega_MELT**2 * (889.190900140*x[0]-698.861205845*self.mu_MELT)

        // dxdt[0] = x[1] / 889.190900140 
        // dxdt[1] = -2 * self.xi_IDLE * self.omega_IDLE * (x[1]) - self.omega_IDLE**2 * (889.190900140 *x[0]-961.996347735*self.mu_IDLE)
        scalar xi_melt = p_const(5);
        scalar omega_melt = p_const(3);
        scalar mu_melt = p_const(1);

        scalar xi_idle = p_const(4);
        scalar omega_idle = p_const(2);
        scalar mu_idle = p_const(0);

        scalar slope = p_const(6);
        scalar offset = p_const(7);

        dxdt(0) =  x(1) / 889.190900140 ;

        dxdt(1) = model_regime * (-2 * xi_melt * omega_melt * (x(1)) - omega_melt*omega_melt * (889.190900140*x(0)-698.861205845*mu_melt) ) + 
                (1-model_regime)*(-2 * xi_idle * omega_idle * (x(1)) - omega_idle*omega_idle * (889.190900140*x(0)-961.996347735*mu_idle) );
        
        dxdt(2) = day_ahead_price * (p_const(8) * 1./(1. + cexp(-(slope * (889.190900140*x(0) - 592.010123492*offset)),15)) +  // Divide by 1000 to get into [0,1]
                                    model_regime * p_const(9) + 
                                    (1-model_regime) * p_const(10)) * 1./60.;                         // Cost
    };
    template <typename scalar>
    scalar Plant::objective(const vector<scalar> &x,
            const vector<scalar> &p_dynamic, const vector<scalar> &p_opt, const vector<double> &p_const) {
        return x(2);
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