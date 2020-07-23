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
  // WALLACE, MISTA
  Curator *wally_{nullptr};
  // Pointer to class that manages GSL functions (RNG, sampling, etc.)
  RandomNumberManagement *gsl_{nullptr};
  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  /* Probability & test statistic trackers*/
  std::map<std::string, double> p_theory_; // Values from SetParameters()
  std::map<std::string, double> p_actual_; // Values stored by event struct
  std::map<std::string, Vec<std::pair<size_t, size_t>>> test_stats_;
  std::map<std::string, Vec<double>> test_ref_;

  /* Auxiliary functions */
  bool equilibrated_{false};
  double scan_window_{10}; // seconds
  double old_density_avg_{0.0};
  double old_density_var_{0.0};
  Vec<double> motor_densities_;
  size_t n_runs_desired_{100}; // std::numeric_limits<size_t>::max()};
  size_t n_runs_recorded_{0};

  size_t step_active_{0};         // KMC step at which motors become active
  bool population_active_{false}; // Switch for motor activity
  bool tethering_active_{false};
  bool lattice_coop_active_{false}; // Switch for lattice deformation effect
  bool lists_up_to_date_{false};    // Whether or not lists are current

  /* Cooperativity mechanisms*/
  int max_neighbs_{2}; // Maximum number of neighbors a single motor can have
  Vec<double> weight_neighbs_bind_;   // Index scheme: [n_neighbs]
  Vec<double> weight_neighbs_unbind_; // Index scheme: [n_neighbs]
  int lattice_cutoff_{0}; // Range of lattice deformation effect (n_sites)
  Vec<double> weight_lattice_bind_;       // Index scheme: [delta] (n_sites)
  Vec<double> weight_lattice_unbind_;     // Index scheme: [delta] (n_sites)
  Vec<double> weight_lattice_bind_max_;   // Index scheme: [n_neighbs]
  Vec<double> weight_lattice_unbind_max_; // Index scheme: [n_neighbs]

  /* KMC event handling */
  Vec<EVENT_T> events_;    // All possible KMC event objects; arbitrary order
  int n_events_to_exe_{0}; // Number of events to execute at each timestep
  Vec<EVENT_T *> events_to_exe_; // List of events to execute at each timestep

  /* Event probabilities; 'avg' refers to poisson events over an interval dt */
  double p_bind_ATP_i_;  // For ATP binding to heads of singly-bound motors
  double p_hydrolyze_;   // For hydrolyzing ATP-bound motor heads
  double p_avg_bind_i_;  // For binding first head of motors from bulk solution
  double p_avg_bind_ii_; // For binding second head of singly-bound motors
  double p_avg_bind_ATP_ii_; // For ATP binding to heads of doubly-bound motors
  double p_avg_unbind_ii_;   // For unbinding rear heads of doubly-bound motors
  double p_avg_unbind_i_;    // For unbinding heads of singly-bound motors

  /* Population size trackers; 'candidates' are targets for poisson events */
  int n_motors_{0}; // Number of motors in system including reservoir; static
  int n_active_{0}; // Number of active (e.g., bound) motors; dynamic
  int n_docked_{0}; // Number of fully-docked heads able to bind
  int n_bound_NULL_i_{0};  // Number of NULL-bound heads of singly-bound motor
  int n_bound_NULL_ii_{0}; // Number of NULL-bound heads of doubly-bound motors
  int n_bound_ATP_{0};     // Number of ATP-bound motor heads
  int n_bound_ADPP_i_{0};  // Number of singly-bound heads able to unbind
  int n_bound_ADPP_ii_{0}; // Number of doubly-bound heads able to unbind

  /* Population trackers -- 'candidates' are targets for possion events */
  Vec<Kinesin> motors_;   // All motors in system including reservoir; static
  Vec<Kinesin *> active_; // Active motors (e.g., bound or tethered); dynamic
  Vec<ENTRY_T> docked_;   // Fully-docked motor heads ready to bind
  Vec<ENTRY_T> bound_NULL_i_;  // NULL-bound heads of singly-bound motors
  Vec<ENTRY_T> bound_NULL_ii_; // NULL-bound heads of doubly-bound motors
  Vec<ENTRY_T> bound_ATP_;     // ATP-bound motor heads
  Vec<ENTRY_T> bound_ADPP_i_;  // ADPP-bound heads of singly-bound motors
  Vec<ENTRY_T> bound_ADPP_ii_; // ADPP-bound heads of doubly-bound motors

private:
  void SetParameters();
  void GenerateMotors();
  void InitializeLists();
  void InitializeTestEnvironment();
  void InitializeEvents();

public:
  KinesinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void ReportProbabilities();

  Kinesin *GetFreeMotor();
  void UpdateLatticeWeights();
  void AddToActive(Kinesin *motor);
  void RemoveFromActive(Kinesin *motor);

  void RunKMC();
  void CheckEquilibration();
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
