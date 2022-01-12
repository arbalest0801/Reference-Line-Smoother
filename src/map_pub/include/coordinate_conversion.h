#ifndef _COORDINATE_CONVERSION_H_
#define _COORDINATE_CONVERSION_H_

#include <array>

class CoordinateConversion{
    public:
    CoordinateConversion() = delete;

    static void cartiesian_to_frenet(const double rs, const double rx, const double ry, const double rtheta, const double rkappa, 
                                     const double rdkappa, const double x,  const double y,  const double v,  const double a, 
                                     const double theta, const double kappa,
                                     std::array<double, 3>* const p_s_condition,
                                     std::array<double, 3>* const P_d_condition);

    static void frenet_to_cartesian(
        const double rs, const double rx, const double ry, const double rtheta, 
        const double rkappa, const double rdkappa,
        std::array<double, 3>& s_cartesian_coord, std::array<double,3>& l_cartesian_coord,
        double* const ptr_x, double* const ptr_y, 
        double* const ptr_theta, double* const ptr_kappa, double* const ptr_v, double* const ptr_a);

};

#endif