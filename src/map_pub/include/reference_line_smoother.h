#pragma once
#include <vector>
#include <utility>
#include "osqp/osqp.h"
#include "Eigen/Eigen"


class ReferenceLineSmoother{
    public:
    ReferenceLineSmoother() = default;
    // ReferenceLineSmoother(double s):s(s){}
    ~ReferenceLineSmoother() = default;
    
    void set_ref_point(const std::vector<std::pair<double, double>> &ref_point){ref_points_ = ref_point;}

    void set_ref_bound_x(const std::vector<double> &ref_bound_x){ref_bound_x_ = ref_bound_x;}

    void set_ref_bound_y(const std::vector<double> &ref_bound_y){ref_bound_y_ = ref_bound_y;}

    void set_weight_s(const double weight_s){weight_smoothness = weight_s;}

    void set_weight_l(const double weight_l){weight_length = weight_l;}

    void set_weight_r(const double weight_r){weight_ref_deviation = weight_r;}

    //求解设置
    void set_max_iter(const int max_iter) { max_iter_ = max_iter; }

    void set_time_limit(const double time_limit) { time_limit_ = time_limit; }

    void set_verbose(const bool verbose) { verbose_ = verbose; }

    void set_scaled_termination(const bool scaled_termination) { scaled_termination_ = scaled_termination;}

    void set_warm_start(const bool warm_start) { warm_start_ = warm_start; }
     
    const std::vector<double>& opt_x(){return x_;}

    const std::vector<double>& opt_y(){return y_;}

    // 求解
    bool Solve();

    private:

    void CalculateKernel(std::vector<c_float>* P_data,    //计算二次项
                         std::vector<c_int>* P_indices,
                         std::vector<c_int>* P_indptr);

    void CalculateOffset(std::vector<c_float>* q);      //计算一次项

    void CalculateAffineConstraint(std::vector<c_float>* A_data,    //计算约束矩阵和上下限
                                std::vector<c_int>* A_indices,
                                std::vector<c_int>* A_indptr);

    bool OptimizeWithOsqp(          //调用osqp求解
    const size_t kernel_dim, const size_t num_affine_constraint,
    std::vector<c_float>* P_data, std::vector<c_int>* P_indices,
    std::vector<c_int>* P_indptr, std::vector<c_float>* A_data,
    std::vector<c_int>* A_indices, std::vector<c_int>* A_indptr,
    std::vector<c_float>* lower_bounds, std::vector<c_float>* upper_bounds,
    std::vector<c_float>* q, std::vector<c_float>* primal_warm_start,
    OSQPData* data, OSQPWorkspace** work, OSQPSettings* settings);

    void SetPrimalWarmStart(std::vector<c_float>* primal_warm_start);

    private:
    //参考点
    std::vector<std::pair<double, double>> ref_points_;
    //偏离参考点误差bound
    std::vector<double> ref_bound_x_;
    std::vector<double> ref_bound_y_;

    double weight_smoothness = 10.0;
    double weight_length = 1.0;
    double weight_ref_deviation = 1.0;

    // Settings of osqp
    int max_iter_ = 4000;
    double time_limit_ = 10.0;
    bool verbose_ = false;
    bool scaled_termination_ = true;
    bool warm_start_ = true;

    //optimize problem size
    int num_of_points_ = 0;
    int num_of_variables_ = 0;
    int num_of_constraints_ = 0;

    //optimized results
    std::vector<double> x_;
    std::vector<double> y_;
};
