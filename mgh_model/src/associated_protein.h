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
  // Indices for lookup tables correspond to x_dist
  Vec<double> cosine_lookup_;
  Vec<double> extension_lookup_;
  // 1st index is n_neighbs (as in PRC1 neighbs); 2nd is x_dist
  Vec<Vec<double>> weight_bind_ii_;
  Vec<Vec<double>> weight_bind_i_teth_;
  Vec<Vec<Vec<double>>> weight_bind_ii_to_teth_;
  Vec<Vec<Vec<double>>> weight_bind_ii_fr_teth_;

  Curator *wally_{nullptr};
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  // Neighbor lists; not to be confused with no. of PRC1 neighbs
  Vec<Tubulin *> neighbor_sites_;
  Vec<Tubulin *> teth_neighbor_sites_;
  Vec<Tubulin *> teth_neighbor_sites_ii_;
  struct Monomer {
    AssociatedProtein *xlink_;
    Tubulin *site_{nullptr};
    std::string state_{std::string("unbound")};
    bool in_scratch_{false};
    void UpdateState();
    int GetPRC1NeighbCount() { return xlink_->GetPRC1NeighbCount(this); }
    int GetDirectionToRest() { return xlink_->GetDirectionToRest(site_); }
    Monomer *GetOtherHead() {
      if (this == &xlink_->head_one_)
        return &xlink_->head_two_;
      else
        return &xlink_->head_one_;
    }
  };

  // see kinesin header for description of variables
  int ID_{-1};
  int speciesID_{1};
  int active_index_{-1};
  int heads_active_{0};
  // For neighbor lists; not to be confused with no. of PRC1 neighbs
  int n_neighbor_sites_{0};
  int n_teth_neighbor_sites_{0};
  int n_teth_neighbor_sites_ii_{0};

  int max_neighbs_{2};

  // x_dist_ is used to index the xlink extensions for lookup
  // e.g. x_dist_ = 0 means an extension of 3 nm (35 - 32)
  int x_dist_{0};      // in no. of sites; can only be 0 or pos.
  int dist_cutoff_{0}; // maximum value x_dist_ can have
  int rest_dist_{0};   // x_dist at which spring extension is ~0

  double r_0_{0.0};
  double k_spring_{0.0};
  double extension_{0.0}; // current extension of xlink in nm
  double cosine_{0.0};    // of angle of xlink w/ respect to MT

  Monomer head_one_{.xlink_ = this}, head_two_{.xlink_ = this};

  bool tethered_{false};
  Kinesin *motor_{nullptr};

private:
  void SetParameters();
  void InitializeLookupTables();
  void InitializeNeighborLists();

public:
  AssociatedProtein();
  void Initialize(system_parameters *parameters, system_properties *properties,
                  int ID);

  Monomer *GetActiveHead();
  Tubulin *GetActiveHeadSite();
  double GetAnchorCoordinate();
  int GetPRC1NeighbCount(Monomer *head);

  void UpdateNeighborSites_II();
  void UpdateNeighborSites_I_Teth();
  void UpdateNeighborSites_II_Teth();
  void UpdateExtension();
  void ForceUnbind();
  void UntetherSatellite();

  int GetDirectionToRest(Tubulin *site);
  double GetExtensionForce(Tubulin *site);
  double GetBindingWeight_II(Tubulin *neighbor);
  double GetBindingWeight_I_Teth(Tubulin *neighbor);
  double GetBindingWeight_II_Teth(Tubulin *neighbor);
  Tubulin *GetSiteCloserToTethRest();
  Tubulin *GetSiteFartherFromTethRest();
  Tubulin *GetWeightedSite_Bind_II();
  Tubulin *GetWeightedSite_Bind_I_Teth();
  Tubulin *GetWeightedSite_Bind_II_Teth();
};
#endif