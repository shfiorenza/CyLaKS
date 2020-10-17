#ifndef _ASSOCIATED_PROTEIN_
#define _ASSOCIATED_PROTEIN_
#include <string>
#include <vector>
class Curator;
class Tubulin;
class Kinesin;
class Microtubule;
struct system_properties;
struct system_parameters;

class AssociatedProtein {
private:
  template <typename DATA_T> using Vec = std::vector<DATA_T>;
  // Neighbor lists for 2nd-stage binding (not to be confused w/ PRC1 neighbs)
  int n_neighbors_bind_ii_{0};
  int n_neighbors_bind_i_teth_{0};
  int n_neighbors_bind_ii_teth_{0};
  Vec<Tubulin *> neighbors_bind_ii_;
  Vec<Tubulin *> neighbors_bind_i_teth_;
  Vec<Tubulin *> neighbors_bind_ii_teth_;
  // Pointers to global system parameters & properties
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};
  // Pointer to system curator, Wallace
  Curator *wally_{nullptr};

public:
  // Monomer structure -- each xlink has two
  struct Monomer {
    AssociatedProtein *xlink_;
    Tubulin *site_{nullptr};
    Monomer *GetOtherHead() {
      if (this == &xlink_->head_one_) {
        return &xlink_->head_two_;
      } else {
        return &xlink_->head_one_;
      }
    }
    int GetDirectionTowardRest() {
      return xlink_->GetDirectionTowardRest(site_);
    }
    int GetPRC1NeighborCount();
  };
  Monomer head_one_{.xlink_ = this}, head_two_{.xlink_ = this};

  // see kinesin header for description of variables
  int id_{-1};
  int species_id_{1};
  int active_index_{-1};
  int heads_active_{0};
  // x_dist_ is used to index the xlink extensions for lookup
  // e.g. x_dist_ = 0 means an extension of 3 nm (35 - 32)
  int x_dist_{0};      // in no. of sites; can only be 0 or pos.
  int dist_cutoff_{0}; // maximum value x_dist_ can have
  int teth_cutoff_{0};
  int comp_cutoff_{0};

  bool tethered_{false};
  Kinesin *motor_{nullptr};

private:
  void SetParameters();
  void InitializeNeighborLists();

public:
  AssociatedProtein();
  void Initialize(system_parameters *, system_properties *, int id);
  // Functions related to baseline function
  Monomer *GetActiveHead();
  double GetAnchorCoordinate();
  void UpdateExtension();
  void ForceUnbind();
  int GetDirectionTowardRest(Tubulin *site);
  // Functions related to tethering
  bool HasSatellite();
  void UntetherSatellite();
  Tubulin *GetSiteCloserToTethRest();
  Tubulin *GetSiteFartherFromTethRest();
  // Functions related to neighbor lists
  void UpdateNeighbors_Bind_II();
  void UpdateNeighbors_Bind_I_Teth();
  void UpdateNeighbors_Bind_II_Teth();
  double GetWeight_Bind_II(Tubulin *neighbor);
  double GetWeight_Bind_I_Teth(Tubulin *neighbor);
  double GetWeight_Bind_II_Teth(Tubulin *neighbor);
  double GetTotalWeight_Bind_II();
  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Bind_II_Teth();
  Tubulin *GetWeightedSite_Bind_II();
  Tubulin *GetWeightedSite_Bind_I_Teth();
  Tubulin *GetWeightedSite_Bind_II_Teth();
};
#endif