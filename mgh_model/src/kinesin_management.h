#ifndef _KINESIN_MANAGEMENT
#define _KINESIN_MANAGEMENT
#include "entry.hpp"
#include "event.hpp"
#include <map>
struct system_parameters;
struct system_properties;
class Curator;
class RandomNumberManagement;

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
  bool lattice_coop_active_{false};
  bool lists_up_to_date_{false};
  // WALLACE, MISTA
  Curator *wally_{nullptr};
  // Pointer to class that manages GSL functions (RNG, sampling, etc.)
  RandomNumberManagement *gsl_{nullptr};
  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  // Index scheme: [tubulin_affinity][n_neighbs][x_dub]
  // If not applicable, will be padded w/ zeros, e.g. [0][0][n_neighbs]
  std::map<std::string, Vec<Vec<Vec<double>>>> p_theory_;
  std::map<std::string, Vec<Vec<Vec<double>>>> p_actual_;
  Vec<std::pair<int, int>> lattice_bind_stats_;
  int test_delta_{-1};
  std::pair<int, int> lattice_step_bind_ii_stats_;
  std::pair<int, int> lattice_step_unbind_ii_stats_;
  std::pair<int, int> lattice_step_unbind_i_stats_;

  // coop stuff -- neighbs
  int max_neighbs_{2};
  // Index scheme: [n_neighbs]
  Vec<double> weight_neighbs_bind_;
  Vec<double> weight_neighbs_unbind_;
  // coop stuff -- lattice
  int lattice_cutoff_{0};
  // Index scheme: [delta], i.e., distance in n_sites
  Vec<double> weight_lattice_bind_;
  Vec<double> weight_lattice_unbind_;
  // Index scheme: [n_neighbs], since each neighb weight has a unique max
  Vec<double> weight_lattice_bind_max_;
  Vec<double> weight_lattice_unbind_max_;

  // Populations are untethered/mixed unless otherwise specified
  int n_motors_{0}; // Total number of motors in system
  int n_active_{0}; // Motors actively bound to some MT/xlink
  // int n_bound_NULL_{0};
  int n_bound_ATP_{0};
  int n_bind_ii_candidates_{0};
  int n_unbind_ii_candidates_{0};
  int n_unbind_i_candidates_{0};
  // Index scheme: [n_neighbs (behind only)]
  Vec<int> n_bound_NULL_;

  // Event probabilities
  double p_hydrolyze_;
  // Avg probabilities used for poisson-based events
  double p_avg_bind_i_;
  double p_avg_bind_ii_;
  double p_avg_unbind_ii_;
  double p_avg_unbind_i_;
  // Index scheme: [n_neighbs (behind only)]
  Vec<double> p_bind_ATP_;

  // 1-D vectors, index is simply motor entry
  Vec<Kinesin> motors_;
  Vec<Kinesin *> active_;
  Vec<ENTRY_T> bound_ATP_;
  // Candidates for poisson-based events
  Vec<ENTRY_T> bind_ii_candidates_;
  Vec<ENTRY_T> unbind_ii_candidates_;
  Vec<ENTRY_T> unbind_i_candidates_;
  // Index scheme: [n_neighbs (behind only)][motor_entry]
  Vec<Vec<ENTRY_T>> bound_NULL_;

private:
  void CalculateCutoffs();
  void SetParameters();
  void GenerateMotors();
  void InitializeLists();
  void InitializeEvents();
  void InitializeTestEnvironment();
  void InitializeTestEvents();

public:
  KinesinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void ReportProbabilities();

  Kinesin *GetFreeMotor();

  void AddToActive(Kinesin *motor);
  void RemoveFromActive(Kinesin *motor);

  void Update_Weights();
  void Update_Docked();
  void Update_Bound_NULL();
  void Update_Bound_ATP();
  void Update_Bound_ADPP_I();
  void Update_Bound_ADPP_II();

  void RunKMC();
  void UpdateLists();
  void SampleEventStatistics();
  void GenerateExecutionSequence();
  void ExecuteEvents();

  void Bind_I(SITE_T *unnoc_site);
  void Bind_ATP(POP_T *bound_head);
  void Hydrolyze(POP_T *bound_head);
  void Bind_II(POP_T *docked_head);
  void Unbind_II(POP_T *bound_head);
  void Unbind_I(POP_T *bound_head);
};
#endif
