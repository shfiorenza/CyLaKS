#ifndef _CYLAKS_ANGULAR_SPRING_HPP_
#define _CYLAKS_ANGULAR_SPRING_HPP_
#include "object.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class AngularSpring : public Object {
private:
  double k_rot_{0.0};

  double torque_{0.0};
  Vec<double> f_vec_;

  Object *rot_point_{nullptr}; // Where torque & force acts
  Object *end_point_{nullptr}; // Other end of rod acting as angular spring

public:
  double theta_min_{0.0};
  double theta_rest_{0.0};
  double theta_max_{0.0};

private:
  void SetCutoffs() {
    // Find cutoff values by setting a max Boltzmann weight of 1e3
    double max_weight{1e3};
    // Recall, Weight = exp(0.5 * E / kbT) [assume lambda = 0.5]
    double E_max{std::log(max_weight) * Params::kbT};
    // E = 0.5 * k_rot * (theta - theta0)^2
    theta_min_ = theta_rest_ - sqrt(2 * E_max / k_rot_);
    theta_max_ = theta_rest_ + sqrt(2 * E_max / k_rot_);
    printf("theta: %g <= %g <= %g\n", theta_min_, theta_rest_, theta_max_);
  }
  void ForceUnbind(Object *obj) { obj->Unbind(); }

public:
  AngularSpring() {}
  void Initialize(size_t sid, size_t id, Object *rotpoint, Object *endpoint,
                  double theta_0, double k_rot) {
    Object::Initialize(sid, id);
    rot_point_ = rotpoint;
    end_point_ = endpoint;
    theta_rest_ = theta_0 * (M_PI / 180.0); // convert to rad
    k_rot_ = k_rot;
    SetCutoffs();
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      f_vec_.push_back(0.0);
    }
  }
  bool UpdatePosition() {
    double r_sq{0.0};
    Vec<double> r_hat(_n_dims_max, 0.0);
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] = rot_point_->pos_[i_dim] - end_point_->pos_[i_dim];
      r_sq += Square(r_hat[i_dim]);
    }
    double r_mag{sqrt(r_sq)};
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] /= r_mag;
    }
    Vec<double> u1{rot_point_->GetBoundObjectOrientation()};
    Vec<double> u2{end_point_->GetBoundObjectOrientation()};
    double theta1{M_PI - acos(Dot(r_hat, u1))};
    // if (theta1 > M_PI / 2) {
    //   theta1 = M_PI - theta1;
    // }
    double theta2{acos(Dot(r_hat, u2))};
    // if (theta2 > M_PI / 2) {
    //   theta2 - M_PI - theta2;
    // }
    printf("r_hat: [%g, %g]\n", r_hat[0], r_hat[1]);
    printf("u1 = [%g, %g]\n", u1[0], u1[1]);
    printf("u2 = [%g, %g]\n", u2[0], u2[1]);
    printf("theta_1/2 = %g / %g\n", theta1, theta2);
    if (theta1 < theta_min_ or theta1 > theta_max_) {
      printf("theta1 = %g!!\n\n", theta1);
      ForceUnbind(rot_point_);
      return false;
    }
    if (theta2 < theta_min_ or theta2 > theta_max_) {
      printf("theta2 = %g!!\n\n", theta2);
      ForceUnbind(end_point_);
      return false;
    }
    double dtheta1{theta1 - theta_rest_};
    double dtheta2{theta2 - theta_rest_};
    double tq1{dtheta1 != 0.0 ? -k_rot_ * (dtheta1 / sin(dtheta1)) : 0.0};
    double tq2{dtheta2 != 0.0 ? -k_rot_ * (dtheta2 / sin(dtheta2)) : 0.0};
    Vec<double> f1_hat{Cross(r_hat, Cross(r_hat, u1))};
    Vec<double> f2_hat{Cross(r_hat, Cross(r_hat, u2))};
    // Torque only comes from dtheta of anchor point
    torque_ = -tq1 * Cross(r_hat, u1);
    if (r_hat[0] > 0.0) {
      torque_ *= -1;
    }
    printf("dtheta_1/2 = %g / %g\n", dtheta1, dtheta2);
    printf("tq1 = %g\n", tq1);
    printf("tq = %g \n", torque_);
    printf("f1_hat: [%g, %g]\n", f1_hat[0], f1_hat[1]);
    printf("f2_hat: [%g, %g]\n\n", f2_hat[0], f2_hat[1]);
    // Force comes from both endpoints due to conservation of angular momentum
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      f_vec_[i_dim] = tq1 / r_mag * f1_hat[i_dim] + tq2 / r_mag * f2_hat[i_dim];
      //   printf("f[i_dim] = %g\n", f_vec_[i_dim]);
    }
    return true;
  }
  void ApplyForces() {
    rot_point_->AddForce(f_vec_);
    rot_point_->AddTorque(torque_);
  }
};

#endif