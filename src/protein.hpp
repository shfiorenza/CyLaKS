#ifndef _CYLAKS_PROTEIN_HPP_
#define _CYLAKS_PROTEIN_HPP_
#include "binding_head.hpp"
#include "linear_spring.hpp"

class Protein : public Object {
protected:
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
  }

  bool HasSatellite();
  void UntetherSatellite();

  void UpdateNeighborList();
  Object *GetWeightedNeighbor();
  BindingHead *GetActiveHead();

  double GetWeight_Bind_I_Teth();
  double GetWeight_Bind_II();
  double GetWeight_Bind_II_Teth();
  BindingSite *GetNeighbor_Bind_I_Teth();
  BindingSite *GetNeighbor_Bind_II();
  BindingSite *GetNeighbor_Bind_II_Teth();

  virtual double GetAnchorCoordinate() {}

  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);
  virtual bool Tether();
  virtual bool Untether();
};

#endif