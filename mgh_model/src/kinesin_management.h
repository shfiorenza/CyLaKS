#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"
#include <string>
#include <functional>

struct system_parameters;
struct system_properties;

//XXX to-do:
//XXX   - make UpdateAllLists() one for loop w/ a shitton of conditions
//				- wait, would this break reproducability??
//					- potentially random order in lists; same RNG
//XXX 	- reorganize KMC functions so non-tether group is first
//			- this will let us do gen_kmc_event in 2 independent pieces
//				- don't need to initialize all arrays if tethers disabled!

class KinesinManagement{
	private:
		// Structure that holds all pertinent info for a given MC event:
		struct event{
			std::string label_;			// Name of event
			std::string target_pop_; 	// Name of pop. this event targets
			int index_ = -1;			// Serialized index of event
			int num_code_ = 0;
			std::function<int(double, int, int)> *sample_stats_;
			int *pop_ptr_ = nullptr; 	// Pointer to pop. size variable
			double p_event_ = 0;		// Probability of event
			int n_events_ = -1;			// Number of predicted events 
			std::function<void(int)> *execute_event_; 
		};
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		int n_motors_ = 0;		// Total number of motors in system
		int n_active_ = 0;		// Motors actively bound to some MT or xlink

		// Populations are untethered/mixed unless otherwise specified 
		int n_free_tethered_ = 0;
		int n_docked_ = 0; 
		int n_bound_NULL_ = 0;
		int n_bound_ATP_ = 0;
		int n_bound_ADPP_i_ = 0;
		int n_bound_ADPP_i_tethered_ = 0; 
		int n_bound_ADPP_ii_ = 0;
		int n_bound_untethered_ = 0;
		std::vector<int> n_bound_tethered_;		// Indexed by x_dub
		
		// See kinesin header for meaningful description of below
		int dist_cutoff_;
		int comp_cutoff_;
		double rest_dist_;

		// Event probabilities 
		double p_bind_i_;
		double p_bind_ATP_;
		double p_hydrolyze_;
		double p_bind_ii_;
		double p_unbind_ii_;
		double p_unbind_i_;
		double p_tether_free_;
		double p_untether_free_;	
		// Below event probabilities are indexed by x_dub
		std::vector<double> p_bind_i_tethered_;
		std::vector<double> p_unbind_i_tethered_;	
		std::vector<double> p_tether_bound_;
		std::vector<double> p_untether_bound_;

		// 1-D vectors, index is simply motor entry
		std::vector<Kinesin> motors_; 
		std::vector<Kinesin*> active_;

		std::vector<Kinesin::head*> docked_; 
		std::vector<Kinesin::head*> bound_NULL_;
		std::vector<Kinesin::head*> bound_ATP_;
		std::vector<Kinesin::head*> bound_ADPP_i_;
		std::vector<Kinesin::head*> bound_ADPP_i_tethered_;
		std::vector<Kinesin::head*> bound_ADPP_ii_;	// FIXME
		std::vector<Kinesin*> free_tethered_;
		std::vector<Kinesin*> bound_untethered_;
		// 2-D vectors, indices are simply [x_dub][motor_entry]
		std::vector< std::vector<Kinesin*> > bound_tethered_;

		std::vector<event> serial_events_;
		std::vector<event> kmc_list_;

	private:
		void GenerateMotors();
		void SetParameters();
		void InitializeLists();
		void InitializeEvents(); 

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		Kinesin* GetFreeMotor();
		Kinesin* GetBoundUntetheredMotor();
		int GetNumBoundUntethered();

		void UpdateAllLists();
		void UpdateFreeTethered();
		void UpdateDocked();
		void UpdateBoundNULL();
		void UpdateBoundATP();
		void UpdateBoundADPP_I();
		void UpdateBoundADPP_I_Tethered();	//XXX
		void UpdateBoundADPP_II();
		void UpdateBoundUntethered();
		void UpdateBoundTethered();

		void GenerateKMCList();
		void UpdateEvents();

		void RunKMC();
		void KMC_Bind_I();	// Bind free ADP head; convert to NULL
		void KMC_Bind_I_Tethered(int x_dub);
		void KMC_Bind_ATP();	// Bind ATP to NULL bound heads
		void KMC_Hydrolyze(); // Convert ATP to ADPP on a bound head
		void KMC_Bind_II();	// Bind docked head; other head must be ADPP
		void KMC_Unbind_II(); 	// Unbind ADPP heads; converts to ADP
		void KMC_Unbind_I();
		void KMC_Unbind_I_Tethered(int x_dub);
		void KMC_Tether_Free();
		void KMC_Tether_Bound(int x_dub);
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dub); 
};
#endif
