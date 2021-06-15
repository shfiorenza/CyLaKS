#ifndef _CYLAKS_LINEAR_SPRING_HPP_
#define _CYLAKS_LINEAR_SPRING_HPP_
#include "object.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class BindingSite;

// LinearSpring: Hookean; works to hold two points at some rest distance
// Note: currently contains some junk related to angular springs; to fix
class LinearSpring : public Object {
protected:
  double k_slack_{0.0};  // For when r < r_0; pN/nm
  double k_spring_{0.0}; // For when r > r_0; pN/nm
  double k_rot_{0.0};

  double dr_{0.0};
  Vec<double> torque_;
  Vec2D<double> f_vec_; // Vector of force for each endpoint (points to center)

  Vec<Object *> endpoints_;

public:
  double r_min_{0.0};  // Minimum spring length (compression); nm
  double r_rest_{0.0}; // Spring rest length; nm
  double r_max_{0.0};  // Maximum spring length (extension); nm

  double theta_min_{0.0};
  double theta_rest_{0.0};
  double theta_max_{0.0};

protected:
  void SetCutoffs() {
    // Find cutoff values by setting a max Boltzmann weight of 1e3
    double max_weight{1e3};
    // Recall, Weight = exp(0.5 * E / kbT) [assume lambda = 0.5]
    double E_max{std::log(max_weight) * Params::kbT};
    // E = 0.5 * k * (r - r0)^2
    r_min_ = r_rest_ - sqrt(2 * E_max / k_slack_);
    r_max_ = r_rest_ + sqrt(2 * E_max / k_spring_);
    // E = 0.5 * k_rot * (theta - theta0)^2
    theta_min_ = theta_rest_ - sqrt(2 * E_max / k_rot_);
    theta_max_ = theta_rest_ + sqrt(2 * E_max / k_rot_);
  }
  void ForceUnbind() {
    // FIXME need to incorporate influence from other springs, e.g. tethers
    if (SysRNG::GetRanProb() < 0.5) {
      endpoints_[0]->Unbind();
    } else {
      endpoints_[1]->Unbind();
    }
  }
  void ForceUnbind(int i_head) { endpoints_[i_head]->Unbind(); }

public:
  LinearSpring() {}
  void Initialize(size_t sid, size_t id, Object *point_one, Object *point_two,
                  double k_slack, double r_0, double k_spring, double theta_0,
                  double k_rot) {
    Object::Initialize(sid, id);
    endpoints_.push_back(point_one);
    endpoints_.push_back(point_two);
    k_slack_ = k_slack;
    r_rest_ = r_0;
    k_spring_ = k_spring;
    theta_rest_ = theta_0 * (M_PI / 180.0); // convert to rad
    k_rot_ = k_rot;
    SetCutoffs();
    torque_.push_back(0.0);
    torque_.push_back(0.0);
    Vec<double> zeroes(_n_dims_max, 0.0);
    f_vec_.push_back(zeroes);
    f_vec_.push_back(zeroes);
  }
  bool UpdatePosition() {
    double r_sq{0.0};                    // Square of spring length
    Vec<double> r_hat(_n_dims_max, 0.0); // points from 2nd to 1st endpoint
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] = endpoints_[0]->pos_[i_dim] - endpoints_[1]->pos_[i_dim];
      r_sq += Square(r_hat[i_dim]);
    }
    double r_mag{sqrt(r_sq)};
    // printf("r = %g\n", r_mag);
    if (r_mag < r_min_ or r_mag > r_max_) {
      // ForceUnbind();
      // return false;
    }
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] /= r_mag;
    }
    // printf("r_hat: [%g, %g]\n", r_hat[0], r_hat[1]);
    // Get theta
    /*
    Vec2D<double> u{endpoints_[0]->GetBoundObjectOrientation(),
                    endpoints_[1]->GetBoundObjectOrientation()};
    Vec<double> theta{M_PI - acos(Dot(r_hat, u[0])), acos(Dot(r_hat, u[1]))};
    // printf("u = [%g, %g] & [%g, %g]\n", u[0][0], u[0][1], u[1][0], u[1][1]);
    // printf("theta = %g & %g\n", theta[0], theta[1]);
    for (int i_endpoint{0}; i_endpoint < 2; i_endpoint++) {
      if (theta[i_endpoint] < theta_min_ or theta[i_endpoint] > theta_max_) {
        // printf("THETA = %g\n", theta[i_endpoint]);
        // ForceUnbind(i_endpoint);
        // return false;
      }
      double dtheta{theta[i_endpoint] - theta_rest_};
      double frac{dtheta != 0.0 ? dtheta / sin(dtheta) : 0.0};
      torque_[i_endpoint] = k_rot_ * frac * Cross(r_hat, u[i_endpoint]);
      // printf("  torque[%i] = %g\n", i_endpoint, torque_[i_endpoint]);
    }
    */
    dr_ = r_mag - r_rest_;
    // printf("dr = %g\n", dr_);
    double f_mag{dr_ > 0.0 ? -k_spring_ * dr_ : -k_slack_ * dr_};
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // Radial forces
      f_vec_[0][i_dim] = f_mag * r_hat[i_dim];
      // Forces from torques
      /*
      f_vec_[0][i_dim] += Cross(r_hat, torque_[0])[i_dim] / r_mag;
      f_vec_[0][i_dim] += Cross(r_hat, torque_[1])[i_dim] / r_mag;
      // printf("  f[0][%i] = %g\n", i_dim, f_vec_[0][i_dim]);
      // Newton's third law; 2nd endpoint gets an equal + opposite force
       */
      f_vec_[1][i_dim] = -f_vec_[0][i_dim];
    }
    return true;
  }
  void ApplyForces() {
    for (int i_endpoint{0}; i_endpoint < endpoints_.size(); i_endpoint++) {
      endpoints_[i_endpoint]->AddForce(f_vec_[i_endpoint]);
      endpoints_[i_endpoint]->AddTorque(torque_[i_endpoint]);
    }
  }
  double GetWeight_Bind(double r) {
    if (r < r_min_ or r > r_max_) {
      return 0.0;
    }
    double dr{r - r_rest_};
    double energy{dr > 0.0 ? 0.5 * k_spring_ * Square(dr)
                           : 0.5 * k_slack_ * Square(dr)};
    return exp(-(1.0 - _lambda_spring) * energy / Params::kbT);
  }
  double GetWeight_Unbind() {
    double energy{dr_ > 0.0 ? 0.5 * k_spring_ * Square(dr_)
                            : 0.5 * k_slack_ * Square(dr_)};
    return exp(_lambda_spring * energy / Params::kbT);
  }
  double GetWeight_Shift(Object *static_site, Object *old_site,
                         Object *new_site) {
    double energy_old{dr_ > 0.0 ? 0.5 * k_spring_ * Square(dr_)
                                : 0.5 * k_slack_ * Square(dr_)};
    double r_x_new{new_site->pos_[0] - static_site->pos_[0]};
    double r_y_new{new_site->pos_[1] - static_site->pos_[1]};
    double r_new{sqrt(Square(r_x_new) + Square(r_y_new))};
    if (r_new < r_min_ or r_new > r_max_) {
      return 0.0;
    }
    double dr_new{r_new - r_rest_};
    double energy_new{dr_new > 0.0 ? 0.5 * k_spring_ * Square(dr_new)
                                   : 0.5 * k_slack_ * Square(dr_new)};
    double dE{energy_new - energy_old};
    // Diffusing towards rest is considered an unbinding-type event in
    // regards to Boltzmann factors, since both events let the spring relax
    if (dE < 0.0) {
      return exp(_lambda_spring * fabs(dE) / Params::kbT);
    }
    // Diffusing away from rest is considered a binding-type event in
    // regards to Boltzmann factors, since both events stretch the spring out
    else {
      return exp(-(1.0 - _lambda_spring) * fabs(dE) / Params::kbT);
    }
  }
};

#endif