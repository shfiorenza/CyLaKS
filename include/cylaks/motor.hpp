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
public:
  CatalyticHead head_one_, head_two_;
  LinearSpring tether_;
  Motor *partner_;

private:
public:
  Motor() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_motor_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_motor_head, this, &head_one_);
    // tether_.Initialize(sid, id, this);
  }
  void ChangeConformation();
  BindingSite *GetDockSite();
  CatalyticHead *GetDockedHead();
  CatalyticHead *GetHeadOne() { return &head_one_; }
  CatalyticHead *GetHeadTwo() { return &head_two_; }
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
  BindingSite *GetNeighbor_Bind_II() { return GetDockSite(); }

  bool UpdateExtension() { return false; }

  void ApplyLatticeDeformation();

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
  bool Tether();
  bool Untether();
};

#endif