#ifndef _CYLAKS_MOTOR_HPP_
#define _CYLAKS_MOTOR_HPP_
#include "catalytic_head.hpp"
#include "linear_spring.hpp"
#include "protein.hpp"

class BindingSite;

class Motor : public Protein {
protected:
  Str ligands_{"yuhh yuh"};

public:
  CatalyticHead head_one_, head_two_;
  LinearSpring tether_;

private:
public:
  Motor(size_t sid, size_t id)
      : Protein(sid, id), head_one_(this, &head_two_, _r_motor_head),
        head_two_(this, &head_one_, _r_motor_head), tether_(this) {}
  BindingSite *GetDockSite();
  CatalyticHead *GetActiveHead();
  void ChangeConformation();
  bool Bind(BindingSite *site, CatalyticHead *head);
  bool Unbind(CatalyticHead *head);
};

#endif