#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_
#define _ASSOCIATED_PROTEIN_MANAGEMENT_
#include "entry.hpp"
#include "event.hpp"
#include <map>
struct system_parameters;
struct system_properties;
class Curator;
class RandomNumberManagement;

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
  // WALLACE, MISTA
  Curator *wally_{nullptr};
  // Pointer to class that manages GSL functions (RNG, sampling, etc.)
  RandomNumberManagement *gsl_{nullptr};
  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  /* Probability & test statistic trackers*/
  // Index scheme: [n_neighbs][x_dub][x]
  // If not applicable, will be padded w/ zeros, e.g. [0][0][n_neighbs]
  std::map<std::string, Vec<Vec<Vec<double>>>> p_theory_;
  std::map<std::string, Vec<Vec<Vec<double>>>> p_actual_;
  std::map<std::string, Vec<std::pair<size_t, size_t>>> test_stats_;
  std::map<std::string, Vec<double>> test_ref_;

  /* System equilibration */
  bool equilibrated_{false};
  double scan_window_{10}; // seconds
  double old_density_avg_{0.0};
  double old_density_var_{0.0};
  Vec<double> densities_;

  /* Population/mechanism activity flags  */
  bool population_active_{false};
  bool crosslinking_active_{false};
  bool tethering_active_{false};
  bool lists_up_to_date_{false};

  /* Energy-dependent effects and their Boltzmann factors */
  int dist_cutoff_{0}; // see assoc. protein header
  int rest_dist_{0};   // see assoc. protein header
  Vec<double> weight_spring_bind_;
  Vec<double> weight_spring_unbind_;
  int max_neighbs_{2};
  Vec<double> weight_neighbs_bind_;   // Index scheme: [n_neighbs]
  Vec<double> weight_neighbs_unbind_; // Index scheme: [n_neighbs]
  int teth_cutoff_{0};
  int comp_cutoff_{0};

  /* KMC event handling */
  Vec<EVENT_T> events_;    // All possible KMC event objects; arbitrary order
  int n_events_to_exe_{0}; // Number of events to execute at each timestep
  Vec<EVENT_T *> events_to_exe_; // List of events to execute at each timestep

  double p_diffuse_i_;
  double p_diffuse_ii_;
  double p_bind_i_;
  double p_bind_i_teth_;
  double p_bind_ii_;
  double p_unbind_ii_;
  double p_unbind_i_;
  double p_tether_;
  double p_untether_;

  /* Population size trackers */
  int n_xlinks_{0}; // Total no. of xlink objects; static
  int n_active_{0}; // No. actively bound; dynamic
  int n_satellites_{0};
  int n_bound_unteth_{0};
  int n_bind_ii_candidates_{0};
  int n_bind_i_teth_candidates_{0};
  int n_bind_ii_teth_candidates_{0};
  Vec<int> n_bound_i_;                      // Indices: [n_neighbs]
  Vec<Vec<int>> n_bound_ii_;                // Indices: [n_neighbs][x]
  Vec<Vec<int>> n_bound_i_teth_;            // Indices: [n_neighbs][x_dub]
  Vec<Vec<Vec<int>>> n_bound_ii_teth_same_; // Indices: [n_neighbs][x_dub][x]
  Vec<Vec<Vec<int>>> n_bound_ii_teth_oppo_; // Indices: [n_neighbs][x_dub][x]

  /* Lists that track different population types */
  Vec<AssociatedProtein> xlinks_; // Actual xlink objects
  Vec<AssociatedProtein *> active_;
  Vec<ENTRY_T> satellites_;
  Vec<ENTRY_T> bound_unteth_;
  Vec<ENTRY_T> bind_ii_candidates_;
  Vec<ENTRY_T> bind_i_teth_candidates_;
  Vec<ENTRY_T> bind_ii_teth_candidates_;
  Vec<Vec<ENTRY_T>> bound_i_;           //   Indices: [n_neighbs][i]
  Vec<Vec<Vec<ENTRY_T>>> bound_ii_;     //   Indices: [n_neighbs][x][i]
  Vec<Vec<Vec<ENTRY_T>>> bound_i_teth_; //   Indices: [n_neighbs][x_dub][i]
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_same_; // [n_neighbs][x_dub][x][i]
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_oppo_; // [n_neighbs][x_dub][x][i]

private:
  void CalculateCutoffs();
  void SetParameters();
  void GenerateXLinks();
  void InitializeLists();
  void InitializeTestEnvironment();
  void InitializeTestEvents();
  void InitializeEvents();

public:
  AssociatedProteinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void ReportProbabilities();

  AssociatedProtein *GetFreeXlink();

  void FlagForUpdate();
  void AddToActive(AssociatedProtein *xlink);
  void RemoveFromActive(AssociatedProtein *xlink);

  void RunKMC();
  void CheckEquilibration();
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
  void Untether(POP_T *satellite_head);
};
#endif