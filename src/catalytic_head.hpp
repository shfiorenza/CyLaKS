#ifndef _CYLAKS_CATALYTIC_HEAD_HPP_
#define _CYLAKS_CATALYTIC_HEAD_HPP_
#include "binding_head.hpp"

enum Ligand { NONE, ATP, ADPP, ADP };

class CatalyticHead : public BindingHead {
private:
  Ligand ligand_{ADP};

  Motor *parent_{nullptr};
  CatalyticHead *other_head_{nullptr};

public:
private:
public:
  CatalyticHead(Object *parent, CatalyticHead *partner, double radius)
      : BindingHead(parent, partner, radius) {}
  virtual ~CatalyticHead();

  Ligand GetLigand() { return ligand_; }
  Motor *GetParent() { return parent_; }
  CatalyticHead *GetOtherHead() { return other_head_; }
  BindingSite *GetSite() { return site_; }
};
#endif