#include <iostream>
#include <limits>
#include "reference_line_smoother.h"
#include "Eigen/Eigen"

using namespace Eigen;
using std::cout;
using std::endl;

template <typename T1, typename T2, int N, typename D>
void DenseToCSCMatrix(const Eigen::Matrix<T1, N, N> &dense_matrix,
                      std::vector<T2> *data, std::vector<D> *indices,
                      std::vector<D> *indptr) {
  static constexpr double epsilon = 1e-9;
  int data_count = 0;
  for (int c = 0; c < dense_matrix.cols(); ++c) {
    indptr->emplace_back(data_count);
    for (int r = 0; r < dense_matrix.rows(); ++r) {
      if (std::fabs(dense_matrix(r, c)) < epsilon) {
        continue;
      }
      data->emplace_back((double)dense_matrix(r, c));
      ++data_count;
      indices->emplace_back(r);
    }
  }
  indptr->emplace_back(data_count);
}

bool ReferenceLineSmoother::Solve(){
    if(ref_points_.empty()) {
        std::cout << "ref points is empty." << std::endl;
        return false;
    }

    if(ref_points_.size() < 3) {
        std::cout << "ref points number smaller than 3." << std::endl;
        return false;
    }

    if(ref_points_.size() != ref_bound_x_.size()) {
        std::cout << "ref point size and bound size not equal." << std::endl;
        return false;
    }

    num_of_points_ = static_cast<int>(ref_points_.size());      
    num_of_variables_ = num_of_points_ * 2;                     //优化变量为(x0,y0,x1,y1,x2,y2......)，同时是kernel矩阵的维度

    // Calculate kernel
    std::vector<c_float> P_data;
    std::vector<c_int> P_indices;
    std::vector<c_int> P_indptr;
    CalculateKernel(&P_data, &P_indices, &P_indptr);        //计算二次项
    // std::cout << "kernel done." << std::endl;

    // Calculate affine constraints
    std::vector<c_float> A_data;
    std::vector<c_int> A_indices;
    std::vector<c_int> A_indptr;

    std::vector<c_float> lower_bounds;
    std::vector<c_float> upper_bounds;
    
    for(int i = 0; i < num_of_points_; ++i){
        upper_bounds.emplace_back(ref_points_[i].first + ref_bound_x_[i]);
        upper_bounds.emplace_back(ref_points_[i].second + ref_bound_y_[i]);
        lower_bounds.emplace_back(ref_points_[i].first - ref_bound_x_[i]);
        lower_bounds.emplace_back(ref_points_[i].second - ref_bound_y_[i]);
    }
    // cout << "bounds done " << endl;
    CalculateAffineConstraint(&A_data, &A_indices, &A_indptr);

    // std::cout << "constraint done." << std::endl;

    // Calculate offset
    std::vector<c_float> q;
    CalculateOffset(&q);        //计算一次项
    // std::cout << "q vector done." << std::endl;

    std::vector<c_float> primal_war_start;

    // SetPrimalWarmStart(&primal_war_start);      //设置原问题warmstart为参考点

    OSQPWorkspace *osqp_workspace = nullptr;        //workspace

    OSQPData *data = reinterpret_cast<OSQPData*>(c_malloc(sizeof(OSQPData)));       //data

    OSQPSettings  *settings = reinterpret_cast<OSQPSettings *>(c_malloc(sizeof(OSQPSettings)));     //settings
    if(settings == nullptr){
        return false;
    } else{
        osqp_set_default_settings(settings);
        settings->max_iter = max_iter_;
        settings->time_limit = time_limit_;
        settings->verbose = verbose_;
        settings->scaled_termination = scaled_termination_;
        settings->warm_start = warm_start_;
        settings->eps_abs = 1.0e-6;
    }
    // cout << "setting done." << endl;

    //前期工作完成，调用求解

    bool result = OptimizeWithOsqp( num_of_variables_, lower_bounds.size(),
                                    &P_data, &P_indices,&P_indptr, 
                                    &A_data, &A_indices, &A_indptr,
                                    &lower_bounds, &upper_bounds,
                                    &q, &primal_war_start,
                                    data, &osqp_workspace, settings);


    if(result == false|| osqp_workspace == nullptr|| osqp_workspace->solution == nullptr){
        // std::cout << "solve failed." << std::endl;
        osqp_cleanup(osqp_workspace);
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
        if(settings) c_free(settings);
        return false;
    }
    //填充xy优化结果
    x_.resize(num_of_points_);  //先resize再填
    y_.resize(num_of_points_);

    for(int i = 0; i < num_of_points_; i++){        //填入求解结果
        int x_index = i * 2;
        x_.at(i) = osqp_workspace->solution->x[x_index];
        y_.at(i) = osqp_workspace->solution->x[x_index + 1];
    }
    // std::cout << "problem solved." << std::endl;
    osqp_cleanup(osqp_workspace);
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
    if(settings) c_free(settings);

    return true;
}

void ReferenceLineSmoother::CalculateKernel(std::vector<c_float>* P_data,    //计算二次项
                         std::vector<c_int>* P_indices,
                         std::vector<c_int>* P_indptr){

                        MatrixXf P1 = MatrixXf::Zero(num_of_variables_, num_of_variables_ - 4);
                        
                        MatrixXf P1_unit = MatrixXf::Zero(6, 2);
                        P1_unit << 1.0, 0.0,
                                   0.0, 1.0,
                                  -2.0, 0.0,
                                   0.0,-2.0,
                                   1.0, 0.0,
                                   0.0, 1.0;     

                        for(int i = 0; i < num_of_points_ - 2; ++i){
                            P1.block(2 * i, 2 * i, 6, 2) = P1_unit;
                        }
                        // cout << "P1 is \n" << P1 << endl;

                        MatrixXf P2 = MatrixXf::Zero(num_of_variables_, num_of_variables_ - 2);
                        MatrixXf P2_unit = MatrixXf::Zero(4, 2);
                        P2_unit << 1.0, 0.0,
                                   0.0, 1.0,
                                  -1.0, 0.0,
                                   0.0,-1.0;

                        for(int i = 0; i < num_of_points_ - 1; ++i){
                            P2.block(2 * i, 2 * i, 4, 2) = P2_unit;
                        }
                        // cout << "P2 is \n" << P2 << endl;

                        MatrixXf P3 = MatrixXf::Identity(num_of_variables_, num_of_variables_);
                        // cout << "p1p1t is \n" << P1 * P1.transpose() << endl;

                        MatrixXf P = weight_smoothness * P1 * P1.transpose() + weight_length * P2 * P2.transpose() + weight_ref_deviation * P3 * P3.transpose();
                        P = P * 2.0;
                        for(int row = 0; row < P.rows(); ++row){
                            for(int col = 0; col < row; ++col){
                                P(row,col) = 0.0;
                            }
                        }
                        // cout << "P is \n" << P << endl;
                        DenseToCSCMatrix(P, P_data, P_indices, P_indptr);
                        // for(const auto& pdata: *P_data){cout << pdata << endl;}
                        // for(const auto& pindices: *P_indices){cout << pindices << endl;}
                        // for(const auto& pind: *P_indptr){cout << pind << endl;}
                        }

void ReferenceLineSmoother::CalculateOffset(std::vector<c_float>* q){//计算一次项
    for(const auto& point : ref_points_){
        q->emplace_back(weight_ref_deviation * -2.0 * point.first);
        q->emplace_back(weight_ref_deviation * -2.0 * point.second);
    }    
}      

void ReferenceLineSmoother::CalculateAffineConstraint(std::vector<c_float>* A_data,    //计算约束矩阵
                                std::vector<c_int>* A_indices,
                                std::vector<c_int>* A_indptr){
                            int ind_A = 0;
                            for(int i = 0; i < num_of_variables_; ++i){
                                A_data->emplace_back(1.0);         //数
                                A_indices->emplace_back(i);        //行号
                                A_indptr->emplace_back(ind_A);    //第几个数
                                ++ind_A;
                            }
                            A_indptr->emplace_back(ind_A);
                        }

void ReferenceLineSmoother::SetPrimalWarmStart(std::vector<c_float>* primal_warm_start){
    for(const auto ref_point : ref_points_){
        primal_warm_start->emplace_back(ref_point.first);
        primal_warm_start->emplace_back(ref_point.second);
    }
}

bool ReferenceLineSmoother::OptimizeWithOsqp(          //调用osqp求解
const size_t kernel_dim, const size_t num_affine_constraint,
std::vector<c_float>* P_data, std::vector<c_int>* P_indices,
std::vector<c_int>* P_indptr, std::vector<c_float>* A_data,
std::vector<c_int>* A_indices, std::vector<c_int>* A_indptr,
std::vector<c_float>* lower_bounds, std::vector<c_float>* upper_bounds,
std::vector<c_float>* q, std::vector<c_float>* primal_warm_start,
OSQPData* data, OSQPWorkspace** work, OSQPSettings* settings){

data->n = kernel_dim;
data->m = num_affine_constraint;

data->P = csc_matrix(data->n, data->n, P_data->size(), P_data->data(), P_indices->data(), P_indptr->data());
data->q = q->data();
data->A = csc_matrix(data->m, data->n, A_data->size(), A_data->data(), A_indices->data(), A_indptr->data());
data->l = lower_bounds->data();
data->u = upper_bounds->data();

// cout << "data done." << endl;

osqp_setup(work, data, settings);

// osqp_warm_start_x(*work, primal_warm_start->data());
// cout << "before solve" << endl;

osqp_solve(*work);

c_int status = (*work)->info->status_val;

if((status != 1 && status != 2) || status < 0){
    // std::cout << "solve failed." << std::endl;
    return false;
}

return true;
}