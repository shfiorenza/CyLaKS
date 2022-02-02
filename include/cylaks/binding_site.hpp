#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"
#include "system_namespace.hpp"

class BindingHead;
class Protofilament;

// BindingSite: Binds various proteins; make up the binding lattice of filaments
class BindingSite : public Sphere {
protected:
  double binding_affinity_{1.0};
  Vec<BindingSite *> neighbors_;
  double weight_bind_{0.0};
  double weight_unbind_{0.0};

public:
  size_t index_{0};
  BindingHead *occupant_{nullptr};
  Protofilament *filament_{nullptr};

private:
public:
  BindingSite() {}
  void Initialize(size_t sid, size_t id, double radius, size_t index,
                  Protofilament *filament) {
    Sphere::Initialize(sid, id, radius);
    index_ = index;
    filament_ = filament;
  }
  void AddNeighbor(BindingSite *site) { neighbors_.emplace_back(site); }
  void SetBindingAffinity(double val) { binding_affinity_ = val; }

  bool IsOccupied() {
    if (occupant_ == nullptr) {
      return false;
    }
    return true;
  }

  void SetWeight_Bind(double val) { weight_bind_ = val; }
  void SetWeight_Unbind(double val) { weight_unbind_ = val; }
  void AddWeight_Bind(double val) {
    int n_neighbs{GetNumNeighborsOccupied()};
    if (weight_bind_ == Sys::weight_lattice_bind_max_[n_neighbs]) {
      return;
    }
    weight_bind_ *= val;
    // If weight is greater than unity, check if it ever gets GREATER THAN max
    if (weight_bind_ > 1.0) {
      if (weight_bind_ > Sys::weight_lattice_bind_max_[n_neighbs]) {
        weight_bind_ = Sys::weight_lattice_bind_max_[n_neighbs];
      }
    }
    // Else if weight is less than unity ,check if it ever gets LESS THAN max
    else if (weight_bind_ < 1.0) {
      if (weight_bind_ < Sys::weight_lattice_bind_max_[n_neighbs]) {
        weight_bind_ = Sys::weight_lattice_bind_max_[n_neighbs];
      }
    }
    // (If weight is equal to unity, neither case matters)
  }
  void AddWeight_Unbind(double val) {
    int n_neighbs{GetNumNeighborsOccupied()};
    if (weight_unbind_ == Sys::weight_lattice_unbind_max_[n_neighbs]) {
      return;
    }
    weight_unbind_ *= val;
    // If weight is greater than unity, check if it ever gets GREATER THAN max
    if (weight_unbind_ > 1.0) {
      if (weight_unbind_ > Sys::weight_lattice_unbind_max_[n_neighbs]) {
        weight_unbind_ = Sys::weight_lattice_unbind_max_[n_neighbs];
      }
    }
    // Else if weight is less than unity ,check if it ever gets LESS THAN max
    else if (weight_unbind_ < 1.0) {
      if (weight_unbind_ < Sys::weight_lattice_unbind_max_[n_neighbs]) {
        weight_unbind_ = Sys::weight_lattice_unbind_max_[n_neighbs];
      }
    }
    // (If weight is equal to unity, neither case matters)
  }
  double GetWeight_Bind() { return weight_bind_ / binding_affinity_; }
  double GetWeight_Unbind() { return weight_unbind_ * binding_affinity_; }

  int GetNumNeighborsOccupied();
  BindingSite *GetNeighbor(int dir);

  void AddForce(Vec<double> f_applied);
  void AddTorque(double tq);
};
#endif