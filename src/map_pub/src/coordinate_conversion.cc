#include "coordinate_conversion.h"
#include <cmath>

void CoordinateConversion::cartiesian_to_frenet(
    const double rs, const double rx, const double ry, const double rtheta, const double rkappa, const double rdkappa,
    const double x,  const double y,  const double v,  const double a, 
    const double theta, const double kappa,
    std::array<double, 3>* const p_s_condition,
    std::array<double, 3>* const p_d_condition){
const double dx = x - rx;
const double dy = y - ry;
p_d_condition->at(0) = 
    std::copysign(std::sqrt(dx * dx + dy * dy), (dy * std::cos(rtheta) - dx * std::sin(rtheta)));

const double one_minus_kl = 1 - rkappa * p_d_condition->at(0);
p_d_condition->at(1) = one_minus_kl * std::tan(theta - rtheta);


const double cos_delta_theta = std::cos(theta - rtheta);
const double kl_minus_kl = rdkappa * p_d_condition->at(0) + rkappa * p_d_condition->at(1);
p_d_condition->at(2) = -1.0 * kl_minus_kl * std::tan(theta - rtheta) + one_minus_kl/(cos_delta_theta/ cos_delta_theta) * 
    ((one_minus_kl)/cos_delta_theta * kappa - rkappa);
    
p_s_condition->at(0) = rs;

p_s_condition->at(1) = v * cos_delta_theta/ one_minus_kl;

const double temp = kappa * (one_minus_kl/cos_delta_theta) - rkappa;

p_s_condition->at(2) = (a * cos_delta_theta - p_s_condition->at(1) * p_s_condition->at(1) * 
    (p_d_condition->at(1) * temp - kl_minus_kl))/ one_minus_kl;
}


void CoordinateConversion::frenet_to_cartesian(
        const double rs, const double rx, const double ry, const double rtheta, 
        const double rkappa, const double rdkappa,
        std::array<double, 3>& s_cartesian_coord, std::array<double,3>& l_cartesian_coord,
        double* const ptr_x, double* const ptr_y, 
        double* const ptr_theta, double* const ptr_kappa, double* const ptr_v, double* const ptr_a){
            
        const double cos_theta_r = std::cos(rtheta);
        const double sin_theta_r = std::sin(rtheta);
        *ptr_x = rx - l_cartesian_coord.at(0) * sin_theta_r;
        *ptr_y = ry + l_cartesian_coord.at(0) * cos_theta_r;


        const double one_minus_kl = 1.0 - rkappa * l_cartesian_coord.at(0);

        const double delta_theta = std::atan2(l_cartesian_coord.at(1), one_minus_kl);

        *ptr_theta = 

        *ptr_v = std::sqrt()

        }

double NormalizeAngle(const double angle){  //-pi to pi



}