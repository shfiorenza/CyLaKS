#ifndef _CYLAKS_LINEAR_SPRING_HPP_
#define _CYLAKS_LINEAR_SPRING_HPP_
#include "object.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class LinearSpring : public Object {
private:
  double k_slack_{0.0};  // For when r < r_0; pN/nm
  double k_spring_{0.0}; // For when r > r_0; pN/nm

  Vec2D<double> f_vec_; // Vector of force for each endpoint (points to center)

  Vec<Object *> endpoints_;

public:
  double r_min_{0.0};  // Minimum spring length (compression); nm
  double r_rest_{0.0}; // Spring rest length; nm
  double r_max_{0.0};  // Maximum spring length (extension); nm

private:
  void SetCutoffs() {
    // Find cutoff values by setting a max Boltzmann weight of 1e3
    double max_weight{1e3};
    // Recall, Weight = exp(0.5 * E / kbT) [assume lambda = 0.5]
    double E_max{std::log(max_weight) * Params::kbT};
    // E = 0.5 * k * (r - r0)^2
    r_min_ = r_rest_ - sqrt(2 * E_max / k_slack_);
    r_max_ = r_rest_ + sqrt(2 * E_max / k_spring_);
  }
  void ForceUnbind() {
    // FIXME need to incorporate influence from other springs, e.g. tethers
    if (SysRNG::GetRanProb() < 0.5) {
      endpoints_[0]->Unbind();
    } else {
      endpoints_[1]->Unbind();
    }
  }

public:
  LinearSpring() {}
  void Initialize(size_t sid, size_t id, Object *point_one, Object *point_two,
                  double k_slack, double r_0, double k_spring) {
    Object::Initialize(sid, id);
    k_slack_ = k_slack;
    r_rest_ = r_0;
    k_spring_ = k_spring;
    SetCutoffs();
    endpoints_.push_back(point_one);
    endpoints_.push_back(point_two);
    Vec<double> zeroes(_n_dims_max, 0.0);
    f_vec_.push_back(zeroes);
    f_vec_.push_back(zeroes);
  }
  bool UpdatePosition() {
    double r_sq{0.0}; // Square of spring length
    // Points from spring center to each endpoint
    Vec2D<double> r_hat(2, Vec<double>(_n_dims_max, 0.0));
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      pos_[i_dim] = Avg(endpoints_[0]->pos_[i_dim], endpoints_[1]->pos_[i_dim]);
      r_sq += Square(endpoints_[0]->pos_[i_dim] - endpoints_[1]->pos_[i_dim]);
      // r_hat vectors are not unit vectors yet; will be rescaled lated
      r_hat[0][i_dim] = endpoints_[0]->pos_[i_dim] - pos_[i_dim];
      r_hat[1][i_dim] = endpoints_[1]->pos_[i_dim] - pos_[i_dim];
    }
    double r_mag{sqrt(r_sq)};
    if (r_mag < r_min_ or r_mag > r_max_) {
      ForceUnbind();
      return false;
    }
    double dr{r_mag - r_rest_};
    double f_mag{dr > 0.0 ? -k_spring_ * dr : -k_slack_ * dr};
    for (int i_endpoint{0}; i_endpoint < endpoints_.size(); i_endpoint++) {
      for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
        // Rescale r_hat so that it is indeed a unit vector
        r_hat[i_endpoint][i_dim] /= (r_mag / 2);
        // Get force vector using f_mag and r_hat
        f_vec_[i_endpoint][i_dim] = f_mag * r_hat[i_endpoint][i_dim];
      }
    }
    return true;
  }
  void ApplyForces() {
    for (int i_endpoint{0}; i_endpoint < endpoints_.size(); i_endpoint++) {
      endpoints_[i_endpoint]->AddForce(f_vec_[i_endpoint]);
    }
  }
};

#endif