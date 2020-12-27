#ifndef _CYLAKS_PROTEIN_HPP_
#define _CYLAKS_PROTEIN_HPP_
#include "binding_head.hpp"
#include "linear_spring.hpp"
#include "system_namespace.hpp"
// FIXME shouldnt include this
#include "binding_site.hpp"
#include "protofilament.hpp"

class Protein : public Object {
protected:
  int dist_cutoff_{10};
  int n_neighbors_bind_ii_{0};
  Vec<BindingSite *> neighbors_bind_ii_;

public:
  size_t active_index_{0};
  size_t n_heads_active_{0};
  BindingHead head_one_, head_two_;
  LinearSpring spring_;

  bool tethered_{false};
  Protein *partner_{nullptr};

private:
  void InitializeNeighborList();
  void ForceUnbind() {}

public:
  Protein() {}
  void Initialize(size_t sid, size_t id) {
    using namespace Params::Xlinks;
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_xlink_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_xlink_head, this, &head_one_);
    spring_.Initialize(sid, id, &head_one_, &head_two_, r_0, k_spring);
    neighbors_bind_ii_.resize(2 * dist_cutoff_ + 1);
  }

  bool HasSatellite();
  void UntetherSatellite();

  int GetNumHeadsActive() { return n_heads_active_; }
  BindingHead *GetActiveHead();
  BindingHead *GetHeadOne() { return &head_one_; }
  BindingHead *GetHeadTwo() { return &head_two_; }

  void UpdateNeighbors_Bind_II();
  double GetWeight_Bind_II(BindingSite *neighb);

  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Bind_II();
  double GetTotalWeight_Bind_II_Teth();
  BindingSite *GetNeighbor_Bind_I_Teth();
  BindingSite *GetNeighbor_Bind_II();
  BindingSite *GetNeighbor_Bind_II_Teth();

  virtual double GetAnchorCoordinate(int i_dim) {
    if (n_heads_active_ != 2) {
      Sys::ErrorExit("Protein::GetAnchorCoord()");
    }
    return (head_one_.site_->pos_[i_dim] + head_two_.site_->pos_[i_dim]) / 2;
  }
  virtual int GetDirectionTowardRest(BindingHead *head) {
    if (n_heads_active_ == 1) {
      return 1;
    } else if (n_heads_active_ == 2) {
      // printf("%g > %g?\n", head->site_->pos_[0], GetAnchorCoordinate(0));
      if (head->site_->pos_[0] > GetAnchorCoordinate(0)) {
        return -1;
      } else if (head->site_->pos_[0] < GetAnchorCoordinate(0)) {
        return 1;
      } else {
        return 0;
      }
    } else {
      Sys::ErrorExit("Protein::GetDirToRest()\n");
    }
    return 0;
  }
  virtual void UpdateExtension() {
    if (n_heads_active_ != 2) {
      return;
    }
    // Update head positions
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      head_one_.pos_[i_dim] = head_one_.site_->pos_[i_dim];
      head_two_.pos_[i_dim] = head_two_.site_->pos_[i_dim];
    }
    // Update spring position
    bool within_cutoff{spring_.UpdatePosition()};
    if (!within_cutoff) {
      ForceUnbind();
      return;
    }
    // If spring is still attached after update, apply forces
    spring_.ApplyForces();
  }
  virtual double GetWeight_Unbind_II(BindingHead *head);
  virtual double GetWeight_Diffuse(BindingHead *head, int dir);

  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);
  virtual bool Diffuse(BindingHead *head, int dir);
  virtual bool Tether();
  virtual bool Untether();
};

#endif