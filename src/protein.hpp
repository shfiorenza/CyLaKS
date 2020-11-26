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
  Object *partner_{nullptr};

private:
  void InitializeNeighborList();

public:
  Protein(size_t sid, size_t id)
      : Object(sid, id), head_one_(this, &head_two_, _r_generic),
        head_two_(this, &head_one_, _r_generic), spring_(this) {}

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

  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);
  virtual bool Tether();
  virtual bool Untether();
};

#endif