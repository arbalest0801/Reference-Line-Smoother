#include "coordinate_conversion.h"
#include <cmath>

void CoordinateConversion::cartiesian_to_frenet(
    const double rs, const double rx, const double ry, const double rtheta, const double rkappa, const double rdkappa,
    const double x,  const double y,  const double v,  const double a, 
    const double theta, const double kappa,
    std::array<double, 3>* const p_s_cartesian_coord,
    std::array<double, 3>* const p_l_cartesian_coord){
const double dx = x - rx;
const double dy = y - ry;
p_l_cartesian_coord->at(0) = 
    std::copysign(std::sqrt(dx * dx + dy * dy), (dy * std::cos(rtheta) - dx * std::sin(rtheta)));

const double one_minus_kl = 1 - rkappa * p_l_cartesian_coord->at(0);
p_l_cartesian_coord->at(1) = one_minus_kl * std::tan(theta - rtheta);


const double cos_delta_theta = std::cos(theta - rtheta);
const double kl_minus_kl = rdkappa * p_l_cartesian_coord->at(0) + rkappa * p_l_cartesian_coord->at(1);
p_l_cartesian_coord->at(2) = -1.0 * kl_minus_kl * std::tan(theta - rtheta) + one_minus_kl/(cos_delta_theta/ cos_delta_theta) * 
    ((one_minus_kl)/cos_delta_theta * kappa - rkappa);
    
p_s_cartesian_coord->at(0) = rs;

p_s_cartesian_coord->at(1) = v * cos_delta_theta/ one_minus_kl;

const double temp = kappa * (one_minus_kl/cos_delta_theta) - rkappa;

p_s_cartesian_coord->at(2) = (a * cos_delta_theta - p_s_cartesian_coord->at(1) * p_s_cartesian_coord->at(1) * 
    (p_l_cartesian_coord->at(1) * temp - kl_minus_kl))/ one_minus_kl;
}


void CoordinateConversion::frenet_to_cartesian(
        const double rs, const double rx, const double ry, const double rtheta, 
        const double rkappa, const double rdkappa,
        std::array<double, 3>& s_cartesian_coord, std::array<double,3>& l_cartesian_coord,
        double* const ptr_x, double* const ptr_y, 
        double* const ptr_theta, double* const ptr_kappa, double* const ptr_v, double* const ptr_a){
            
        const double cos_theta_r = std::cos(rtheta);
        const double sin_theta_r = std::sin(rtheta);

        *ptr_x = rx - sin_theta_r * l_cartesian_coord[0];
        *ptr_y = ry + cos_theta_r * l_cartesian_coord[0];

        const double one_minus_kappa_r_d = 1 - rkappa * l_cartesian_coord[0];

        const double tan_delta_theta = l_cartesian_coord[1] / one_minus_kappa_r_d;
        const double delta_theta = std::atan2(l_cartesian_coord[1], one_minus_kappa_r_d);
        const double cos_delta_theta = std::cos(delta_theta);

        *ptr_theta = NormalizeAngle(delta_theta + rtheta);

        const double kappa_r_d_prime =
            rdkappa * l_cartesian_coord[0] + rkappa * l_cartesian_coord[1];
        *ptr_kappa = (((l_cartesian_coord[2] + kappa_r_d_prime * tan_delta_theta) *
                        cos_delta_theta * cos_delta_theta) /
                            (one_minus_kappa_r_d) +
                        rkappa) *
                    cos_delta_theta / (one_minus_kappa_r_d);

        const double d_dot = l_cartesian_coord[1] * s_cartesian_coord[1];
        *ptr_v = std::sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d *
                                s_cartesian_coord[1] * s_cartesian_coord[1] +
                            d_dot * d_dot);

        const double delta_theta_prime =
            one_minus_kappa_r_d / cos_delta_theta * (*ptr_kappa) - rkappa;

        *ptr_a = s_cartesian_coord[2] * one_minus_kappa_r_d / cos_delta_theta +
                s_cartesian_coord[1] * s_cartesian_coord[1] / cos_delta_theta *
                    (l_cartesian_coord[1] * delta_theta_prime - kappa_r_d_prime);
}

double NormalizeAngle(const double angle){  //-pi to pi
    double a = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (a < 0.0) 
    {
        a += (2.0 * M_PI);
    }
    return a - M_PI;
}