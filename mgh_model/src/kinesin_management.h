#ifndef _KINESIN_MANAGEMENT
#define _KINESIN_MANAGEMENT
#include "entry.hpp"
#include "event.hpp"
#include <map>
struct system_parameters;
struct system_properties;
class Curator;

class KinesinManagement {
private:
  // Use a template for 'Vec' rather than std::vector (aesthetic only)
  template <class DATA_T> using Vec = std::vector<DATA_T>;
  // Define data types
  using POP_T = Kinesin::Monomer;
  using ALT_T = AssociatedProtein::Monomer;
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
  double t_active_{0.0};
  bool population_active_{false};
  bool tethering_active_{false};
  bool lists_up_to_date_{false};
  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};
  // WALLACE, MISTA
  Curator *wally_{nullptr};

public:
  // Index scheme: [tubulin_affinity][n_neighbs][x_dub]
  // If not applicable, will be padded w/ zeros, e.g. [0][0][n_neighbs]
  std::map<std::string, Vec<Vec<Vec<double>>>> p_theory_;
  std::map<std::string, Vec<Vec<Vec<double>>>> p_actual_;

  // coop stuff
  int n_affinities_{10};
  int n_stacks_{5};
  int n_affinities_tot_{n_affinities_ * n_stacks_ + 1};
  int max_neighbs_{2};

  // See kinesin header for meaningful description of below
  int teth_cutoff_;
  int comp_cutoff_;
  double rest_dist_;

  // Populations are untethered/mixed unless otherwise specified
  int n_motors_{0}; // Total number of motors in system
  int n_active_{0}; // Motors actively bound to some MT/xlink
  int n_bound_ATP_{0};
  int n_bound_ATP_st_{0};
  int n_bound_unteth_{0};
  int n_satellites_{0};
  // Candidates for poisson-based events
  int n_bind_i_teth_candidates_{0};
  int n_tether_bound_candidates_{0};
  // Index scheme: [x_dub]
  Vec<int> n_bound_teth_;
  // Index scheme: [n_neighbs]
  Vec<int> n_bound_NULL_;
  // Index scheme: [n_neighbs (behind only)][x_dub]
  Vec<Vec<int>> n_bound_NULL_to_teth_;
  Vec<Vec<int>> n_bound_NULL_fr_teth_;
  // Index scheme: [tubulin_affinity][n_neighbs]
  Vec<Vec<int>> n_docked_;
  Vec<Vec<int>> n_bound_ADPP_ii_;
  Vec<Vec<int>> n_bound_ADPP_i_;
  Vec<Vec<int>> n_bound_ADPP_i_st_;
  // Index scheme: [tubulin_affinity][n_neighbs][x_dub]
  Vec<Vec<Vec<int>>> n_bound_ADPP_i_teth_;
  Vec<Vec<Vec<int>>> n_bound_ADPP_i_teth_st_;

  // Event probabilities
  double p_hydrolyze_;
  double p_hydrolyze_st_;
  double p_tether_free_;
  double p_untether_satellite_;
  // Avg probabilities used for poisson-based events
  double p_avg_bind_i_teth_;
  double p_avg_tether_bound_;
  // Index scheme: [x_dub]
  Vec<double> p_untether_bound_;    // curent x_dub
  Vec<double> weight_tether_bound_; // proposed x_dub
  // Index scheme: [n_neighbs (behind only)]
  Vec<double> p_bind_ATP_;
  // Index scheme: [n_neighbs (behind only)][x_dub]
  Vec<Vec<double>> p_bind_ATP_to_teth_;
  Vec<Vec<double>> p_bind_ATP_fr_teth_;
  // Index scheme: [tubulin_affinity][n_neighbs]
  Vec<Vec<double>> p_bind_i_;
  Vec<Vec<double>> p_bind_ii_;
  Vec<Vec<double>> p_unbind_ii_;
  Vec<Vec<double>> p_unbind_i_;
  Vec<Vec<double>> p_unbind_i_st_;
  // Index scheme: [tubulin_affinity][n_neighbs][x_dub]
  Vec<Vec<Vec<double>>> p_unbind_i_teth_;    // current x_dub
  Vec<Vec<Vec<double>>> p_unbind_i_teth_st_; // current x_dub
  Vec<Vec<Vec<double>>> weight_bind_i_teth_; // proposed x_dub

  // 1-D vectors, index is simply motor entry
  Vec<Kinesin> motors_;
  Vec<Kinesin *> active_;
  Vec<ENTRY_T> bound_ATP_;
  Vec<ENTRY_T> bound_ATP_st_;
  Vec<ENTRY_T> bound_unteth_;
  Vec<ENTRY_T> satellites_;
  // Candidates for poisson-based events
  Vec<ENTRY_T> bind_i_teth_candidates_;
  Vec<ENTRY_T> tether_bound_candidates_;
  // Index scheme: [x_dub][motor_entry]
  Vec<Vec<ENTRY_T>> bound_teth_;
  // Index scheme: [n_neighbs (behind only)][motor_entry]
  Vec<Vec<ENTRY_T>> bound_NULL_;
  // Index scheme: [n_neighbs (behind only)][x_dub][motor_entry]
  Vec<Vec<Vec<ENTRY_T>>> bound_NULL_to_teth_;
  Vec<Vec<Vec<ENTRY_T>>> bound_NULL_fr_teth_;
  // Index scheme: [tubulin_affinity][n_neighbs][motor_entry]
  Vec<Vec<Vec<ENTRY_T>>> docked_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_ii_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_i_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_i_st_;
  // Index scheme: [tubulin_affinity][n_neighbs][x_dub][motor_entry]
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ADPP_i_teth_;
  Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ADPP_i_teth_st_;

private:
  void CalculateCutoffs();
  void SetParameters();
  void GenerateMotors();
  void InitializeLists();
  void InitializeEvents();

public:
  KinesinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void ReportProbabilities();

  Kinesin *GetFreeMotor();

  void AddToActive(Kinesin *motor);
  void RemoveFromActive(Kinesin *motor);

  void ReportExecutionOf(std::string event_name);
  void ReportFailureOf(std::string event_name);

  void FlagForUpdate();

  void Update_Extensions();
  void Update_Docked();
  void Update_Bound_NULL();
  void Update_Bound_NULL_Teth();
  void Update_Bound_ATP();
  void Update_Bound_ATP_Stalled();
  void Update_Bound_ADPP_I();
  void Update_Bound_ADPP_I_Stalled();
  void Update_Bound_ADPP_I_Teth();
  void Update_Bound_ADPP_I_Teth_Stalled();
  void Update_Bound_ADPP_II();
  void Update_Bound_Unteth();
  void Update_Bound_Teth();
  void Update_Free_Teth();

  double GetWeight_Bind_I_Teth();
  double GetWeight_Tether_Bound();

  void Set_Bind_I_Teth_Candidates(int n_to_set);
  void Set_Tether_Bound_Candidates(int n_to_set);

  void RunKMC();
  void UpdateLists();
  void SampleEventStatistics();
  void GenerateExecutionSequence();
  void ExecuteEvents();

  void Bind_I(SITE_T *unnoc_site);
  void Bind_I_Teth(POP_T *satellite_head);
  void Bind_ATP(POP_T *bound_head);
  void Hydrolyze(POP_T *bound_head);
  void Bind_II(POP_T *docked_head);
  void Unbind_II(POP_T *bound_head);
  void Unbind_I(POP_T *bound_head);
  void Tether_Free(ALT_T *untethered_head);
  void Tether_Bound(POP_T *bound_head);
  void Untether(POP_T *head);
};
#endif
