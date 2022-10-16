#ifndef _CYLAKS_BINDING_HEAD_HPP_
#define _CYLAKS_BINDING_HEAD_HPP_
#include "sphere.hpp"

class BindingSite;
class Protein;

// BindingHead: Can bind to and diffuse along binding sites
class BindingHead : public Sphere {
protected:
public:
  BindingSite *site_{nullptr};

  Protein *parent_{nullptr};
  BindingHead *other_head_{nullptr};

private:
public:
  BindingHead() {}
  void Initialize(size_t sid, size_t id, double radius, Protein *parent_ptr,
                  BindingHead *other_head_ptr) {
    Sphere::Initialize(sid, id, radius);
    parent_ = parent_ptr;
    other_head_ = other_head_ptr;
  }
  void Initialize(size_t sid, size_t id, double radius) {
    Sphere::Initialize(sid, id, radius);
  }

  virtual bool Trailing() { return false; }

  int GetDirectionTowardRest();
  virtual int GetNumHeadsActive();
  virtual int GetNumNeighborsOccupied();
  virtual int GetNumNeighborsOccupied_Side();

  BindingHead *GetOtherHead() { return other_head_; }
  BindingSite *GetSite() { return site_; }

  virtual Vec<double> GetSpringOrientation() {
    double r_sq{0.0};
    Vec<double> r_hat(_n_dims_max, 0.0);
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] = pos_[i_dim] - GetOtherHead()->pos_[i_dim];
      r_sq += Square(r_hat[i_dim]);
    }
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      r_hat[i_dim] /= sqrt(r_sq);
    }
    return r_hat;
  }
  virtual Vec<double> GetBoundObjectOrientation();

  virtual void AddForce(Vec<double> f_applied);
  virtual void AddTorque(double tq);

  virtual bool UntetherSatellite();

  virtual double GetWeight_Unbind_II();
  virtual bool Unbind();

  virtual double GetWeight_Diffuse(int dir);
  virtual bool Diffuse(int dir);
  virtual bool Diffuse_Side(int dir);
};
#endif