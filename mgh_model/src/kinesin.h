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
  // Indices for lookup tables correspond to distance
  // in (no. of sites)/2, NOT extension of tether
  // e.g. ...lookup_[1] is an x-dist of 1/2 of a site
  // Index scheme: [x_dist_doubled] (current)
  Vec<double> cosine_lookup_;
  Vec<double> extension_lookup_;
  // Index scheme: [x_dub] (proposed)
  Vec<double> weight_tether_bound_;
  // Index scheme: [tubulin_affinity][n_kif4a_neighbs][x_dub] (proposed)
  Vec<Vec<Vec<double>>> weight_bind_i_teth_;
  // Neighbor lists for tethering (not to be confused w/ Kif4a neighbs)
  int n_neighbors_bind_i_teth_{0};
  int n_neighbors_tether_bound_{0};
  Vec<Tubulin *> neighbors_bind_i_teth_;
  Vec<AssociatedProtein *> neighbors_tether_bound_;
  // Pointers to global system parameters & properties
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};
  // Pointer to system curator, Wallace
  Curator *wally_{nullptr};

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
    void RelieveFrustration() {
      if (motor_->frustrated_) {
        if (trailing_) {
          motor_->ChangeConformation();
        } else {
          motor_->frustrated_ = false;
        }
      }
    }
    int GetCometSize() { return motor_->GetCometSize(this); }
    int GetAffinity();
    int GetKif4ANeighborCount();
    int GetKif4ANeighborCount_Step();
  };

  int id_;               // Unique id of this kinesin in resevoir
  int species_id_{2};    // Unique id of this species (kinesin)
  int active_index_{-1}; // index of this motor in active_ list
  int heads_active_{0};
  // x_dist_dub is used to index the tether extension of motors, e.g.
  // an x_dist_dub of 10 means an extension of -80 nm (40 - 120)
  int x_dist_doubled_{0}; // in no. of sites
  int comp_cutoff_{0};    // min value x_dist (not 2x) can be
  double rest_dist_{0};   // spring extension is ~0 for this
  int teth_cutoff_{0};    // max value x_dist (not 2x) can be

  double r_0_{0.0};       // in nm
  double k_spring_{0.0};  // in pN / nm
  double k_slack_{0.0};   // in pN / nm
  double cosine_{0.0};    // of motor tether angle w.r.t. horizontal
  double extension_{0.0}; // in nm

  Monomer head_one_{.motor_ = this}, head_two_{.motor_ = this};

  bool frustrated_{false};
  bool tethered_{false};
  Microtubule *mt_{nullptr};
  AssociatedProtein *xlink_{nullptr};

private:
  void SetParameters();
  void InitializeLookupTables();
  void InitializeNeighborLists();

public:
  Kinesin();
  void Initialize(system_parameters *, system_properties *, int id);
  // Functions related to baseline function
  Monomer *GetActiveHead();
  Monomer *GetDockedHead();
  Tubulin *GetDockSite();
  double GetDockedCoordinate();
  int GetCometSize(Monomer *head);
  void ChangeConformation();
  bool IsStalled();
  // Functions related to tethering
  int GetDirectionTowardRest();
  double GetStalkCoordinate();
  double GetRestLengthCoordinate();
  double GetTetherForce(Tubulin *site);
  void UpdateExtension();
  void ForceUntether();
  void UntetherSatellite();
  bool HasSatellite(); // FIXME
  bool AtCutoff();
  // Functions related to neighbor lists
  void UpdateNeighbors_Bind_I_Teth();
  void UpdateNeighbors_Tether_Bound();
  double GetWeight_Bind_I_Teth(Tubulin *site);
  double GetWeight_Bind_II();
  double GetWeight_Unbind_II();
  double GetWeight_Tether_Bound(AssociatedProtein *xlink);
  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Tether_Bound();
  Tubulin *GetWeightedSite_Bind_I_Teth();
  AssociatedProtein *GetWeightedXlink_Tether_Bound();
};
#endif