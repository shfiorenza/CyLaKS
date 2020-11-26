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
  // Neighbor lists for tethering (not to be confused w/ Kif4a neighbs)
  int n_neighbors_bind_i_teth_{0};
  int n_neighbors_tether_bound_{0};
  Vec<Tubulin *> neighbors_bind_i_teth_;
  Vec<AssociatedProtein *> neighbors_tether_bound_;
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
  // x_dist_dub is used to index the tether extension of motors, e.g.
  // an x_dist_dub of 10 means an extension of -80 nm (40 - 120)
  int x_dist_doubled_{0}; // in no. of sites
  int teth_cutoff_{0};
  int comp_cutoff_{0};

  bool tethered_{false};
  Microtubule *mt_{nullptr};
  AssociatedProtein *xlink_{nullptr};

private:
  void SetParameters();
  void InitializeNeighborLists();

public:
  Kinesin();
  void Initialize(system_parameters *, system_properties *, int id);
  // Functions related to baseline function
  void ChangeConformation();
  double GetWeight_Bind_II();
  double GetWeight_Unbind_II();
  double GetWeight_BindATP_II();
  Monomer *GetActiveHead();
  Monomer *GetDockedHead();
  Tubulin *GetDockSite();
  // Functions related to tethering
  void UpdateExtension(); // FIXME
  void ForceUntether();
  void UntetherSatellite();
  bool HasSatellite();
  bool IsStalled();
  bool IsSteppingTowardRest();
  int GetDirectionTowardRest();
  double GetStalkCoordinate();
  double GetRestLengthCoordinate();
  double GetTetherForce(Tubulin *site);
  // Functions related to neighbor lists
  void UpdateNeighbors_Bind_I_Teth();
  void UpdateNeighbors_Tether_Bound();
  double GetWeight_Bind_I_Teth(Tubulin *site);
  double GetWeight_Tether_Bound(AssociatedProtein *xlink);
  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Tether_Bound();
  Tubulin *GetWeightedSite_Bind_I_Teth();
  AssociatedProtein *GetWeightedXlink_Tether_Bound();
};
#endif