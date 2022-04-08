#ifndef _CYLAKS_PROTEIN_HPP_
#define _CYLAKS_PROTEIN_HPP_
#include "angular_spring.hpp"
#include "binding_head.hpp"
#include "linear_spring.hpp"
#include "system_namespace.hpp"

class Motor;

// Protein: Essentially just a passive crosslinker at this stage
//          Has two binding heads connected by a linear hookean spring
class Protein : public Object {
protected:
  double ran_{0.0};
  int n_neighbors_bind_ii_{0};
  Vec<BindingSite *> neighbors_bind_ii_;
  int n_neighbors_bind_i_teth_{0};
  Vec<BindingSite *> neighbors_bind_i_teth_;
  int n_neighbors_bind_ii_teth_{0};
  Vec<BindingSite *> neighbors_bind_ii_teth_;

public:
  size_t active_index_{0};
  size_t n_heads_active_{0};

  BindingHead head_one_, head_two_;
  LinearSpring spring_;

  Protein *teth_partner_{nullptr};

protected:
  void InitializeNeighborLists();

  void UpdateNeighbors_Bind_II();
  void UpdateNeighbors_Bind_I_Teth();
  void UpdateNeighbors_Bind_II_Teth();
  double GetSoloWeight_Bind_II(BindingSite *neighb);
  double GetSoloWeight_Bind_I_Teth(BindingSite *target);
  double GetSoloWeight_Bind_II_Teth(BindingSite *neighb);

public:
  Protein() {}
  void Initialize(size_t sid, size_t id) {
    using namespace Params;
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_xlink_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_xlink_head, this, &head_one_);
    spring_.Initialize(sid, id, &head_one_, &head_two_, Xlinks::k_spring,
                       Xlinks::r_0, Xlinks::k_spring, Xlinks::theta_0,
                       Xlinks::k_rot);
    // Maximum possible x_distance of spring will occur when r_y = 0
    size_t x_max{(size_t)std::ceil(spring_.r_max_ / Filaments::site_size)};
    neighbors_bind_ii_.resize(2 * x_max + 1);
    neighbors_bind_i_teth_.resize(Filaments::count *
                                  (2 * Sys::teth_x_max_ + 1));
  }

  Protein *GetTethPartner() { return teth_partner_; }

  int GetNumHeadsActive() { return n_heads_active_; }
  virtual BindingHead *GetActiveHead() {
    if (n_heads_active_ != 1) {
      Sys::ErrorExit("Protein::GetActiveHead()");
    }
    if (head_one_.site_ != nullptr) {
      return &head_one_;
    } else if (head_two_.site_ != nullptr) {
      return &head_two_;
    } else {
      Sys::ErrorExit("AssociatedProtein::GetActiveHead");
      return nullptr;
    }
  }
  virtual BindingHead *GetHeadOne() { return &head_one_; }
  virtual BindingHead *GetHeadTwo() { return &head_two_; }

  bool IsTethered() {
    if (teth_partner_ != nullptr) {
      return true;
    }
    return false;
  }
  bool HasSatellite();
  bool UntetherSatellite();
  double GetWeight_Bind_I_Teth();
  BindingSite *GetNeighbor_Bind_I_Teth();

  virtual void ApplyLatticeDeformation() {}

  virtual bool UpdateExtension();
  virtual int GetDirectionTowardRest(BindingHead *head);
  virtual double GetAnchorCoordinate(int i_dim);

  virtual BindingSite *GetNeighbor_Bind_II();

  virtual double GetWeight_Diffuse(BindingHead *head, int dir);
  virtual double GetWeight_Bind_II();
  virtual double GetWeight_Unbind_II(BindingHead *head);

  virtual bool Diffuse(BindingHead *head, int dir);
  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);
  virtual bool Tether(Protein *teth_partner);
  virtual bool Untether();
};

#endif