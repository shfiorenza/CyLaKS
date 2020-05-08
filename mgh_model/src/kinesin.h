#ifndef _KINESIN
#define _KINESIN
#include <string>
#include <vector>
class Curator;
class Tubulin;
class AssociatedProtein;
class Microtubule;
struct system_properties;
struct system_parameters;

class Kinesin {
private:
  template <typename DATA_T> using Vec = std::vector<DATA_T>;
  // Pointer to system curator, Wallace
  Curator *wally_{nullptr};
  // Pointers to global system parameters & properties
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  // Monomer structure -- each motor has two
  struct Monomer {
    Kinesin *motor_;
    Tubulin *site_{nullptr};
    bool trailing_;
    std::string ligand_{"ADP"};
    Monomer *GetOtherHead() {
      if (this == &motor_->head_one_) {
        return &motor_->head_two_;
      } else {
        return &motor_->head_one_;
      }
    }
    int GetKif4ANeighborCount();
  };
  Monomer head_one_{.motor_ = this}, head_two_{.motor_ = this};

  int id_;               // Unique id of this kinesin in resevoir
  int species_id_{2};    // Unique id of this species (kinesin)
  int active_index_{-1}; // index of this motor in active_ list
  int heads_active_{0};

  Microtubule *mt_{nullptr};

private:
public:
  Kinesin();
  void Initialize(system_parameters *, system_properties *, int id);
  // Functions related to baseline function
  void ChangeConformation();
  double GetWeight_Bind_II();
  double GetWeight_Unbind_II();
  Monomer *GetActiveHead();
  Monomer *GetDockedHead();
  Tubulin *GetDockSite();
};
#endif