#ifndef _KINESIN
#define _KINESIN
#include <string>
#include <vector>
class Tubulin;
class Microtubule;
class AssociatedProtein;
struct system_properties;
struct system_parameters;

class Kinesin {
private:
  // Indices for lookup tables correspond to distance
  // in (no. of sites)/2, NOT extension of tether
  // e.g. ...lookup_[1] is an x-dist of 1/2 of a site
  std::vector<double> cosine_lookup_;
  std::vector<double> extension_lookup_;
  std::vector<double> weight_lookup_;
  std::vector<double> weight_alt_lookup_;
  std::vector<double> p_step_to_rest_;
  std::vector<double> p_step_fr_rest_;

  // Neighbor sites are for when the motor is tethered but unbound
  std::vector<Tubulin *> neighbor_sites_;
  // Neighbor xlinks are for when the motor is bound but untethered
  std::vector<AssociatedProtein *> neighbor_xlinks_;

  system_parameters *parameters_ = nullptr;
  system_properties *properties_ = nullptr;

public:
  struct head {
    Kinesin *motor_;
    Tubulin *site_;
    bool trailing_;
    bool in_scratch_{false};
    std::string ligand_;
    std::string state_;
    Tubulin *stored_dock_site_{nullptr};

    head(Kinesin *parent, Tubulin *site, bool trailing, std::string ligand,
         std::string state)
        : motor_{parent}, site_{site}, trailing_{trailing}, ligand_{ligand},
          state_{state} {}
    head *GetOtherHead() {
      if (this == &motor_->head_one_)
        return &motor_->head_two_;
      else
        return &motor_->head_one_;
    }
  };

  int id_;            // Unique id of this kinesin in resevoir
  int speciesID_ = 2; // Unique id of this species (kinesin)
  int active_index_;  // index of this motor in active_ list
  int heads_active_ = 0;
  int n_neighbor_sites_ = 0;
  int n_neighbor_xlinks_ = 0;

  bool frustrated_ = false;
  bool tethered_ = false;

  // x_dist_dub is used to index the tether extension of motors, e.g.
  // an x_dist_dub of 10 means an extension of -80 nm (40 - 120)
  int x_dist_doubled_; // in no. of sites
  int dist_cutoff_;    // max value x_dist (not 2x) can be
  int comp_cutoff_;    // min value x_dist (not 2x) can be
  double rest_dist_;   // spring extension is ~0 for this

  double r_0_;       // in nm
  double k_spring_;  // in pN / nm
  double k_slack_;   // in pN / nm
  double extension_; // in nm
  double cosine_;    // of motor tether angle w.r.t. horizontal

  // stats from perpetually applied force
  double applied_p_step_;

  head head_one_ = head(this, nullptr, false, "ADP", "unbound"),
       head_two_ = head(this, nullptr, true, "ADP", "unbound");

  AssociatedProtein *xlink_ = nullptr;
  Microtubule *mt_ = nullptr;

private:
  void SetParameters();
  void CalculateCutoffs();
  void InitializeNeighborLists();
  void InitializeWeightLookup();
  void InitializeExtensionLookup();
  void InitializeSteppingProbabilities();

public:
  Kinesin();
  void Initialize(system_parameters *parameters, system_properties *properties,
                  int id);

  head *GetActiveHead();
  head *GetDockedHead();
  head *StoreDockSite();
  double GetStalkCoordinate(); // tail originates from stalk
  double GetDockedCoordinate();
  void ChangeConformation(); // FIXME add applied force w/ teth
  bool IsStalled();

  bool AtCutoff();
  void UpdateNeighborSites();
  void UpdateNeighborXlinks();
  void UpdateExtension();
  void ForceUntether();
  void UntetherSatellite();

  int GetDirectionTowardRest();
  double GetRestLengthCoordinate(); // coord where ext ~ 0 when bound
  double GetTetherForce(Tubulin *site);
  double GetTotalBindingWeight();
  double GetBindingWeight(Tubulin *site);
  double GetTotalTetheringWeight();
  double GetTetheringWeight(AssociatedProtein *xlink);
  Tubulin *GetWeightedNeighborSite();
  AssociatedProtein *GetWeightedNeighborXlink();
};
#endif