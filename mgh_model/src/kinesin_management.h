#ifndef _KINESIN_MANAGEMENT
#define _KINESIN_MANAGEMENT
#include "entry.h"
#include "event.hpp"
class Curator;
struct system_parameters;
struct system_properties;

class KinesinManagement {
private:
  using POP_T = Kinesin::head;
  using ALT_T = AssociatedProtein::Monomer;
  using SITE_T = Tubulin;
  using MGMT_T = KinesinManagement;
  // ENTRY_T is defined in entry.h header
  using EVENT_T = Event<MGMT_T *, ENTRY_T>;
  // Use a template for 'Vec' rather than std::vector (aesthetic)
  template <class DATA_T> using Vec = std::vector<DATA_T>;
  using Str = std::string;
  // All possible KMC event objects; arbitrary sequential order
  Vec<EVENT_T> events_;
  // Evets segregated by target pop. (for stat correction)
  Vec<Vec<EVENT_T *>> events_by_pop_;
  // List of event to execute any given timestep; dynamic
  Vec<EVENT_T *> events_to_exe_;
  // Temporarily holds entries after KMC events
  int n_scratched_ = 0;
  Vec<POP_T *> scratch_;

  bool lists_up_to_date_{false};

  // Pointers to global system params & props; same for all classes
  system_parameters *parameters_ = nullptr;
  system_properties *properties_ = nullptr;
  // WALLACE, MISTA
  Curator *wally_;

public:
  bool verbose_{false};

  // coop stuff
  int n_affinities_{11};
  // int deformation_range_{};
  int max_neighbs_{2};

  // See kinesin header for meaningful description of below
  int teth_cutoff_;
  int comp_cutoff_;
  double rest_dist_;

  // Populations are untethered/mixed unless otherwise specified
  int n_motors_ = 0; // Total number of motors in system
  int n_active_ = 0; // Motors actively bound to some MT/xlink
  int n_free_tethered_ = 0;
  int n_bound_NULL_ = 0;
  int n_bound_ATP_ = 0;
  int n_bound_ATP_stalled_ = 0;
  int n_bound_untethered_ = 0;
  // Below population sizes are indexed by [tubulin_affinity][n_neighbs]
  Vec<Vec<int>> n_docked_;
  Vec<Vec<int>> n_bound_ADPP_i_;
  Vec<Vec<int>> n_bound_ADPP_i_stalled_;
  Vec<Vec<int>> n_bound_ADPP_ii_;
  // Below population sizes are indexed by x_dub
  Vec<int> n_docked_tethered_;
  Vec<int> n_bound_NULL_tethered_;
  Vec<int> n_bound_ADPP_i_tethered_;
  Vec<int> n_bound_ADPP_i_teth_st_;
  Vec<int> n_bound_tethered_;

  // Event probabilities
  double p_bind_i_tethered_;
  double p_bind_ATP_;
  double p_hydrolyze_;
  double p_hydrolyze_stalled_;
  double p_tether_free_;
  double p_tether_bound_;
  double p_untether_free_;
  // Below event probabilities are indexed by [tubulin_affinity][n_neighbs]
  Vec<Vec<double>> p_bind_i_;
  Vec<Vec<double>> p_bind_ii_;
  Vec<Vec<double>> p_unbind_ii_;
  Vec<Vec<double>> p_unbind_i_;
  Vec<Vec<double>> p_unbind_i_stalled_;
  // Below event probabilities are indexed by x_dub
  Vec<double> p_bind_ATP_tethered_;
  Vec<double> p_bind_ii_tethered_;
  Vec<double> p_unbind_i_tethered_;
  Vec<double> p_unbind_i_teth_st_;
  Vec<double> p_untether_bound_;

  // 1-D vectors, index is simply motor entry
  Vec<Kinesin> motors_;
  Vec<Kinesin *> active_;
  Vec<ENTRY_T> free_tethered_;
  Vec<ENTRY_T> bound_untethered_;
  Vec<ENTRY_T> bound_NULL_;
  Vec<ENTRY_T> bound_ATP_;
  Vec<ENTRY_T> bound_ATP_stalled_;
  // 3-D vectors, inidices are [tubulin_affinity][n_neighbs][motor_entry]
  Vec<Vec<Vec<ENTRY_T>>> docked_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_i_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_i_stalled_;
  Vec<Vec<Vec<ENTRY_T>>> bound_ADPP_ii_;
  // 2-D vectors, indices are simply [x_dub][motor_entry]
  Vec<Vec<ENTRY_T>> docked_tethered_;
  Vec<Vec<ENTRY_T>> bound_NULL_tethered_;
  Vec<Vec<ENTRY_T>> bound_ADPP_i_tethered_;
  Vec<Vec<ENTRY_T>> bound_ADPP_i_teth_st_;
  Vec<Vec<ENTRY_T>> bound_tethered_;

private:
  void GenerateMotors();
  void SetParameters();
  void InitializeLists();
  void InitializeEvents();

public:
  KinesinManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);

  Kinesin *GetFreeMotor();

  void Update_All_Lists();
  void Update_Free_Teth();
  void Update_Docked();
  void Update_Docked_Teth();
  void Update_Bound_Teth();
  void Update_Bound_Unteth();
  void Update_Bound_NULL();
  void Update_Bound_NULL_Teth();
  void Update_Bound_ATP();
  void Update_Bound_ATP_Stalled();
  void Update_Bound_ADPP_I();
  void Update_Bound_ADPP_I_Stalled();
  void Update_Bound_ADPP_I_Tethered();
  void Update_Bound_ADPP_I_Tethered_Stalled();
  void Update_Bound_ADPP_II();

  void Run_KMC();
  void Refresh_Populations();
  void Generate_Execution_Sequence();
  int Sample_Event_Statistics();

  double GetWeight_TetherBound();
  double GetWeight_BindTethered();

  POP_T *CheckScratchFor(std::string pop);
  void SaveToScratch(POP_T *head);
  void ReportExecutionOf(std::string event_name);
  void ReportFailureOf(std::string event_name);

  // void Execution_Relay(ENTRY_T target, int code);
  void KMC_Bind_I(SITE_T *unnoc_site);   // Bind free ADP head; convert to NULL
  void KMC_Bind_ATP(POP_T *bound_head);  // Bind ATP to NULL-bound heads
  void KMC_Hydrolyze(POP_T *bound_head); // Convert ATP to ADPP on a bound head
  void KMC_Bind_II(POP_T *docked_head);  // Bind docked head
  void KMC_Unbind_II(POP_T *bound_head); // Unbind ADPP heads; converts to ADP
  void KMC_Unbind_I(POP_T *bound_head);
  //  void KMC_Bind_ATP_Tethered(int x_dub);
  //  void KMC_Bind_I_Tethered();
  // void KMC_Bind_II_Tethered(int x_dub);
  //  void KMC_Unbind_I_Tethered(int x_dub);
  // void KMC_Unbind_I_Tethered_Stalled(int x_dub); // XXX
  //  void KMC_Tether_Free();
  //  void KMC_Tether_Bound();
  //  void KMC_Untether_Free();
  //  void KMC_Untether_Bound(int x_dub);
};
#endif
