#ifndef _CYLAKS_LINEAR_SPRING_HPP_
#define _CYLAKS_LINEAR_SPRING_HPP_
#include "object.hpp"

class LinearSpring : public Object {
private:
  double f_{0.0};
  double r_{0.0};
  Vec<double> r_vec_;

  double r_min_{0.0};  // Minimum spring length (compression); nm
  double r_rest_{0.0}; // Spring rest length; nm
  double r_max_{0.0};  // Maximum spring length (extension); nm

  double k_slack_{0.0};  // For when r < r_0; pN/nm
  double k_spring_{0.0}; // For when r > r_0; pN/nm

  Vec<Object *> endpoints_;

public:
private:
  void ForceUnbind();

public:
  LinearSpring() {}
  void Initialize(size_t sid, size_t id, Object *point_one, Object *point_two,
                  double r_0, double k_spring) {
    Object::Initialize(sid, id);
    r_rest_ = r_0;
    k_spring_ = k_slack_ = k_spring;
    endpoints_.push_back(point_one);
    endpoints_.push_back(point_two);
    r_vec_.resize(_n_dims_max);
  }
  bool UpdatePosition() {
    double r_sq{0.0};
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_vec_[i_dim] = endpoints_[0]->pos_[i_dim] - endpoints_[1]->pos_[i_dim];
      pos_[i_dim] = Avg(endpoints_[0]->pos_[i_dim], endpoints_[1]->pos_[i_dim]);
      r_sq += Square(r_vec_[i_dim]);
    }
    r_ = sqrt(r_sq);
    // if (r_ < r_min_ or r_ > r_max_) {
    //   return false;
    // }
    double dr{r_ - r_rest_};
    if (dr > 0.0) {
      f_ = -k_spring_ * dr;
    } else {
      f_ = -k_slack_ * dr;
    }
    return true;
  }
  void ApplyForces() {
    for (int i_endpoint{0}; i_endpoint < endpoints_.size(); i_endpoint++) {
      // Since r_ points from the 2nd to 1st endpoint, it will be in the same
      // direction as spring displacement for the 1st head. For the 2nd head,
      // we must correct with a factor of -1, essentially flipping the vector
      double correction_factor{i_endpoint == 0 ? 1 : -1};
      Vec<double> f_applied(_n_dims_max, 0.0);
      for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
        f_applied[i_dim] = f_ * correction_factor * (r_vec_[i_dim] / r_);
      }
      endpoints_[i_endpoint]->AddForce(f_applied);
    }
  }
};

#endif