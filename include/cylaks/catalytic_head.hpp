#ifndef _CYLAKS_CATALYTIC_HEAD_HPP_
#define _CYLAKS_CATALYTIC_HEAD_HPP_
#include "binding_head.hpp"

class BindingSite;
class Motor;

// CatalyticHead: A binding head that also undergoes a ATP hydrolysis cycle
//                Cycles thru ADP, NONE, ATP, and ADPP ligands in that order
class CatalyticHead : public BindingHead {
private:
public:
  // enum Ligand { NONE, ATP, ADPP, ADP };
  Ligand ligand_{ADP};

  bool trailing_{false};

  Motor *parent_{nullptr};
  CatalyticHead *other_head_{nullptr};

  CatalyticHead *test_partner_{nullptr}; // ! for lattice_step test mode

private:
public:
  CatalyticHead() {}
  void Initialize(size_t sid, size_t id, double radius, Motor *parent_ptr,
                  CatalyticHead *other_head_ptr);

  int GetNumNeighborsOccupied();
  int GetNumHeadsActive();
  bool Trailing() { return trailing_; }
  Ligand GetLigand() { return ligand_; }
  CatalyticHead *GetOtherHead() { return other_head_; }

  double GetWeight_Unbind_II();
  bool Unbind();
  double GetWeight_Diffuse(int dir);
  bool Diffuse(int dir);
};
#endif