#ifndef _CYLAKS_PROTEIN_HPP_
#define _CYLAKS_PROTEIN_HPP_
#include "binding_head.hpp"
#include "linear_spring.hpp"

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

public:
  Protein() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_xlink_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_xlink_head, this, &head_one_);
    spring_.Initialize(sid, id, this);
    neighbors_bind_ii_.resize(2 * dist_cutoff_ + 1);
  }

  bool HasSatellite();
  void UntetherSatellite();

  int GetNumHeadsActive() { return n_heads_active_; }
  BindingHead *GetActiveHead();
  BindingHead *GetHeadOne() { return &head_one_; }

  // void UpdateNeighborList();
  // Object *GetWeightedNeighbor();

  void UpdateNeighbors_Bind_II();
  double GetWeight_Bind_II(BindingSite *neighb);

  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Bind_II();
  double GetTotalWeight_Bind_II_Teth();
  BindingSite *GetNeighbor_Bind_I_Teth();
  BindingSite *GetNeighbor_Bind_II();
  BindingSite *GetNeighbor_Bind_II_Teth();

  virtual double GetAnchorCoordinate() {}

  virtual void UpdateExtension() {}

  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);
  virtual bool Diffuse(BindingHead *head, int dir);

  virtual bool Tether();
  virtual bool Untether();
};

#endif