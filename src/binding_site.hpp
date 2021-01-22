#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"
#include "system_namespace.hpp"

class BindingHead;
class Protofilament;

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

  void SetBindingAffinity(double val) { binding_affinity_ = val; }

  void AddNeighbor(BindingSite *site) { neighbors_.emplace_back(site); }

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
    // printf("**heyyyy\n");
    if (weight_bind_ == Sys::weight_lattice_bind_max_[n_neighbs]) {
      return;
    }
    // printf("**HEY - val = %g\n", val);
    weight_bind_ *= val;
    // printf("**AH\n");
    if (weight_bind_ > Sys::weight_lattice_bind_max_[n_neighbs]) {
      weight_bind_ = Sys::weight_lattice_bind_max_[n_neighbs];
      // printf("**u wot\n");
      return;
    }
  }
  void AddWeight_Unbind(double val) {
    int n_neighbs{GetNumNeighborsOccupied()};
    // printf("heyyyy\n");
    if (weight_unbind_ == Sys::weight_lattice_unbind_max_[n_neighbs]) {
      return;
    }
    // printf("HEY - val = %g\n", val);
    weight_unbind_ *= val;
    // printf("AH\n");
    if (weight_unbind_ > Sys::weight_lattice_unbind_max_[n_neighbs]) {
      weight_unbind_ = Sys::weight_lattice_unbind_max_[n_neighbs];
      // printf("u wot\n");
      return;
    }
  }
  double GetWeight_Bind() { return weight_bind_ / binding_affinity_; }
  double GetWeight_Unbind() { return weight_unbind_ * binding_affinity_; }

  int GetNumNeighborsOccupied() {
    if (_n_neighbs_max == 0) {
      return 0;
    }
    int n_neighbs{0};
    for (auto const &site : neighbors_) {
      // printf("hi\n");
      if (site->occupant_ != nullptr) {
        n_neighbs++;
      }
    }
    // printf("noh - %i\n", n_neighbs);
    return n_neighbs;
  }

  void AddForce(Vec<double> f_applied);
  void AddTorque(double tq);

  BindingSite *GetNeighbor(int dir);
};
#endif