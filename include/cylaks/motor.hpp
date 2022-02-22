#ifndef _CYLAKS_MOTOR_HPP_
#define _CYLAKS_MOTOR_HPP_
#include "catalytic_head.hpp"
#include "linear_spring.hpp"
#include "protein.hpp"

class BindingSite;

// Motor: Steps via coordinated mechanochemcial cyles of two catalytic heads
//        Can exert forces and drag cargo via its 'tether' (or stalk)
class Motor : public Protein {
protected:
  int n_neighbors_tether_{0};
  Vec<Protein *> neighbors_tether_;

public:
  CatalyticHead head_one_, head_two_;
  LinearSpring tether_;

protected:
  void UpdateNeighbors_Tether();
  double GetSoloWeight_Tether(Protein *xlink);

public:
  Motor() {}
  void Initialize(size_t sid, size_t id) {
    using namespace Params;
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_motor_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_motor_head, this, &head_one_);
    // ! FIXME -- how to use endpoints for tether?
    // ! FIXME -- create Anchor obj. that dynamically updates
    tether_.Initialize(sid, id, &head_one_, &head_two_, Motors::k_slack,
                       Motors::r_0, Motors::k_tether, 0.0, 0.0);
    size_t x_max{(size_t)std::ceil(tether_.r_max_ / Filaments::site_size)};
    neighbors_tether_.resize(Filaments::count * (2 * x_max + 1));
    neighbors_bind_i_teth_.resize(Filaments::count * (2 * x_max + 1));
  }

  CatalyticHead *GetActiveHead() {
    if (n_heads_active_ != 1) {
      Sys::ErrorExit("Motor::GetActiveHead [0]");
      return nullptr;
    }
    if (head_one_.site_ != nullptr) {
      return &head_one_;
    } else if (head_two_.site_ != nullptr) {
      return &head_two_;
    } else {
      Sys::ErrorExit("Motor::GetActiveHead [1]");
      return nullptr;
    }
  }
  CatalyticHead *GetHeadOne() { return &head_one_; }
  CatalyticHead *GetHeadTwo() { return &head_two_; }
  CatalyticHead *GetDockedHead();
  BindingSite *GetDockSite();

  void ChangeConformation();

  BindingSite *GetNeighbor_Bind_II() { return GetDockSite(); }

  void ApplyLatticeDeformation();

  // ! FIXME!
  bool UpdateExtension() { return false; }
  int GetDirectionTowardsRest(CatalyticHead *head);
  void ForceUntether();

  double GetWeight_Tether();
  Protein *GetNeighbor_Tether();

  double GetAnchorCoordinate(int i_dim);

  double GetWeight_Diffuse(CatalyticHead *head, int dir);
  double GetWeight_Bind_II();
  double GetWeight_BindATP_II(CatalyticHead *head);
  double GetWeight_Unbind_II(CatalyticHead *head);
  double GetWeight_Unbind_I();

  bool Diffuse(CatalyticHead *head, int dir);
  bool Bind(BindingSite *site, CatalyticHead *head);
  bool Bind_ATP(CatalyticHead *head);
  bool Hydrolyze(CatalyticHead *head);
  bool Unbind(CatalyticHead *head);
  // bool Tether(Protein *teth_partner);
  // bool Untether();
};

#endif