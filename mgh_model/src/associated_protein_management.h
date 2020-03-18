#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_
#define _ASSOCIATED_PROTEIN_MANAGEMENT_
#include "entry.hpp"
#include "event.hpp"
#include <map>
struct system_parameters;
struct system_properties;
class Curator;

class AssociatedProteinManagement {
private:
  // Use a template for 'Vec' rather than std::vector (aesthetic only)
  template <class DATA_T> using Vec = std::vector<DATA_T>;
  // Define data types
  using POP_T = AssociatedProtein::Monomer;
  using ALT_T = Kinesin::Monomer;
  using SITE_T = Tubulin;
  // ENTRY_T is defined in entry.h header
  using EVENT_T = Event<ENTRY_T>;
  // All possible KMC event objects; arbitrary sequential order
  Vec<EVENT_T> events_;
  // Number of events to execute at any given timestep
  int n_events_to_exe_{0};
  // List of events to execute any given timestep; dynamically updated
  Vec<EVENT_T *> events_to_exe_;
  int verbosity_{0};
  bool population_active_{false};
  bool crosslinking_active_{false};
  bool tethering_active_{false};
  bool lists_up_to_date_{false};
  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};
  // WALLACE, MISTA
  Curator *wally_{nullptr};

public:
  // Index scheme: [n_neighbs][x_dub][x]
  // If not applicable, will be padded w/ zeros, e.g. [0][0][n_neighbs]
  std::map<std::string, Vec<Vec<Vec<double>>>> p_theory_;
  std::map<std::string, Vec<Vec<Vec<double>>>> p_actual_;
  Vec<std::pair<unsigned long, unsigned long>> bind_ii_stats_;
  Vec<std::pair<unsigned long, unsigned long>> bind_i_teth_stats_;
  Vec<Vec<std::pair<unsigned long, unsigned long>>> bind_ii_to_teth_stats_;
  Vec<Vec<std::pair<unsigned long, unsigned long>>> bind_ii_fr_teth_stats_;

  // Neighbor coop stuff; still kinda preliminary
  int max_neighbs_{2};
  double interaction_energy_{0.0}; // in kBT

  int dist_cutoff_{0}; // see assoc. protein header
  int rest_dist_{0};   // see assoc. protein header
  int teth_cutoff_{0}; // see kinesin header (dist_cutoff_ there)
  int comp_cutoff_{0}; // see kinesin header (comp_cutoff_ there)

  /* Population size trackers */
  int n_xlinks_{0}; // Total no. of xlink objects; static
  int n_active_{0}; // No. actively bound; dynamic
  int n_free_teth_{0};
  int n_bound_unteth_{0};
  int n_bind_ii_candidates_{0};
  int n_bind_i_teth_candidates_{0};
  int n_bind_ii_teth_candidates_{0};
  // First index is number of PRC1 neighbors: [0], [1], or [2]
  Vec<int> n_bound_i_;
  // Second index is [x] or [x_dub] for base or teth pops.
  Vec<Vec<int>> n_bound_i_teth_;
  Vec<Vec<int>> n_bound_ii_;
  // Second index is [x_dub], third & final index is [x]
  Vec<Vec<Vec<int>>> n_bound_ii_teth_same_;
  Vec<Vec<Vec<int>>> n_bound_ii_teth_oppo_;

  /* Probabilities of possible KMC events */
  double p_avg_bind_i_teth_;
  double p_avg_bind_ii_;
  double p_tether_free_;
  double p_untether_free_;
  // First index is number of PRC1 neighbors: [0], [1], or [2]
  Vec<double> p_bind_i_;
  Vec<double> p_unbind_i_;
  Vec<double> p_diffuse_i_fwd_;
  Vec<double> p_diffuse_i_bck_;
  // Indices are [n_neighbs][x] or [n_neighbs][x_dub]
  Vec<Vec<double>> p_unbind_i_teth_;
  Vec<Vec<double>> p_unbind_ii_;
  Vec<Vec<double>> p_diffuse_i_to_teth_rest_;
  Vec<Vec<double>> p_diffuse_i_fr_teth_rest_;
  Vec<Vec<double>> p_diffuse_ii_to_rest_;
  Vec<Vec<double>> p_diffuse_ii_fr_rest_;
  // Indices are [n_neighbs][x_dub][x]
  Vec<Vec<Vec<double>>> p_unbind_ii_to_teth_;
  Vec<Vec<Vec<double>>> p_unbind_ii_fr_teth_;
  Vec<Vec<Vec<double>>> p_diffuse_ii_to_both_;
  Vec<Vec<Vec<double>>> p_diffuse_ii_fr_both_;
  Vec<Vec<Vec<double>>> p_diffuse_ii_to_self_fr_teth_;
  Vec<Vec<Vec<double>>> p_diffuse_ii_fr_self_to_teth_;

  Vec<Vec<double>> weight_bind_ii_;     // [n_neighbs][x]
  Vec<Vec<double>> weight_bind_i_teth_; // [n_neighbs][x_dub]
  // x is PROPOSED x_dist after bind_ii; x_dub is current value before bind_ii
  Vec<Vec<Vec<double>>> weight_bind_ii_to_teth_; // [n_neighbs][x_dub][x]
  Vec<Vec<Vec<double>>> weight_bind_ii_fr_teth_; // [n_neighbs][x_dub][x]

  /* Lists that track different population types */
  Vec<AssociatedProtein> xlinks_; // Actual xlink objects
  Vec<AssociatedProtein *> active_;
  Vec<ENTRY_T> free_teth_;
  Vec<ENTRY_T> bound_unteth_;
  // Target pools for poisson processes (do not track n_neighbs/etc.)
  Vec<ENTRY_T> bind_ii_candidates_;
  Vec<ENTRY_T> bind_i_teth_candidates_;
  Vec<ENTRY_T> bind_ii_teth_candidates_;
  // First index is number of PRC1 neighbors: [0], [1], or [2]
  // Second index is actual xlink entry
  Vec<Vec<ENTRY_T>> bound_i_;
  // Second index is [x] or [x_dub]; third index is xlink entry
  Vec<Vec<Vec<ENTRY_T>>> bound_ii_;
  Vec<Vec<Vec<ENTRY_T>>> bound_i_teth_;
  // Second index is [x_dub]; third index is [x]; fourth is entry
  // e.g., [0][16][2][1] -> 2nd xlink w/ x_dub=16, x=2, & 0 neighbs
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_oppo_;
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_same_;

private:
  void CalculateCutoffs();
  void SetParameters();
  void GenerateXLinks();
  void InitializeLists();
  void InitializeEvents();
  void InitializeTestEnvironment();
  void InitializeTestEvents();

public:
  AssociatedProteinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void ReportProbabilities();

  AssociatedProtein *GetFreeXlink();

  void AddToActive(AssociatedProtein *xlink);
  void RemoveFromActive(AssociatedProtein *xlink);

  void FlagForUpdate();

  void Update_Extensions();
  void Update_Bound_I();
  void Update_Bound_I_Teth();
  void Update_Bound_II();
  void Update_Bound_II_Teth();
  void Update_Bound_Unteth();
  void Update_Free_Teth();

  void RunKMC();
  void UpdateLists();
  void SampleEventStatistics();
  void GenerateExecutionSequence();
  void ExecuteEvents();

  void Diffuse(POP_T *head, int dir);
  void Bind_I(SITE_T *target_site);
  void Bind_I_Teth(POP_T *satellite_head);
  void Bind_II(POP_T *bound_head);
  void Unbind_II(POP_T *bound_head);
  void Unbind_I(POP_T *bound_head);
  void Tether_Free(ALT_T *untethered_head);
  void Untether_Free(POP_T *satellite_head);
};
#endif