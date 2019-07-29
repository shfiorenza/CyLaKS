#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "entry.h"
#include <string>
#include <functional>
class Curator; 
struct system_parameters;
struct system_properties;

class KinesinManagement{
	private:
		using POP_T = Kinesin::head;
		Curator* wally_;
		// Structure that holds all pertinent info for a given MC event:
		struct event{
			event(int i, int code, std::string l, std::string t_p,
					std::function<int(double, int)> p_dist, 
					int *pop, double p): index_(i), kmc_code_(code),
					label_(l), target_pop_(t_p), prob_dist_(p_dist), 
					pop_ptr_(pop), p_occur_(p) {}
			void SampleStatistics(){
				if(*pop_ptr_ > 0)
					n_expected_ = prob_dist_(p_occur_, *pop_ptr_); 
				else n_expected_ = 0; 
			}
		private:
		public:
			std::function<int(double, int)> prob_dist_;
			int index_ = -1;			// Serialized unique index of event
			int kmc_code_ = -1;			// Encodes which KMC_ funct to call
			std::string label_ = "BRUH";			// Name of event
			std::string target_pop_; 	// Name of pop. this event targets
			double p_occur_ = 0;		// Probability of event occuring
			int *pop_ptr_ = nullptr; 	// Pointer to pop. size variable
			int n_expected_ = 0;		// Number of predicted events 
		};
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		// Populations are untethered/mixed unless otherwise specified 
		int n_motors_ = 0;		// Total number of motors in system
		int n_active_ = 0;		// Motors actively bound to some MT/xlink
		int n_free_tethered_ = 0;
		int n_docked_ = 0; 
		int n_bound_NULL_ = 0;
		int n_bound_ATP_ = 0;
		int n_bound_ATP_stalled_ = 0; 
		int n_bound_ADPP_i_ = 0;
		int n_bound_ADPP_i_stalled_ = 0; 
		int n_bound_ADPP_ii_ = 0;
		int n_bound_untethered_ = 0;
		// Below population sizes are indexed by x_dub
		std::vector<int> n_docked_tethered_;
		std::vector<int> n_bound_NULL_tethered_;
		std::vector<int> n_bound_ADPP_i_tethered_; 
		std::vector<int> n_bound_ADPP_i_teth_st_;
		std::vector<int> n_bound_tethered_;
		
		// See kinesin header for meaningful description of below
		int dist_cutoff_;
		int comp_cutoff_;
		double rest_dist_;

		// Event probabilities 
		double p_bind_i_;
		double p_bind_i_tethered_;
		double p_bind_ATP_;
		double p_hydrolyze_;
		double p_hydrolyze_stalled_;
		double p_bind_ii_;
		double p_unbind_ii_;
		double p_unbind_i_;
		double p_unbind_i_stalled_;
		double p_tether_free_;
		double p_tether_bound_;
		double p_untether_free_;	
		// Below event probabilities are indexed by x_dub
		std::vector<double> p_bind_ATP_tethered_; 
		std::vector<double> p_bind_ii_tethered_;
		std::vector<double> p_unbind_i_tethered_;	
		std::vector<double> p_unbind_i_teth_st_;
		std::vector<double> p_untether_bound_;

		// 1-D vectors, index is simply motor entry
		std::vector<Kinesin> motors_; 
		std::vector<Kinesin*> active_;
		std::vector<Kinesin*> free_tethered_;
		std::vector<ENTRY_T> bound_untethered_;
		std::vector<Kinesin::head*> docked_; 
		std::vector<Kinesin::head*> bound_NULL_;
		std::vector<Kinesin::head*> bound_ATP_;
		std::vector<Kinesin::head*> bound_ATP_stalled_;
		std::vector<Kinesin::head*> bound_ADPP_i_;
		std::vector<Kinesin::head*> bound_ADPP_i_stalled_;
		std::vector<Kinesin::head*> bound_ADPP_ii_;
		// 2-D vectors, indices are simply [x_dub][motor_entry]
		std::vector< std::vector<Kinesin::head*> > docked_tethered_;
		std::vector< std::vector<Kinesin::head*> > bound_NULL_tethered_;
		std::vector< std::vector<Kinesin::head*> > bound_ADPP_i_tethered_;
		std::vector< std::vector<Kinesin::head*> > bound_ADPP_i_teth_st_;
		std::vector< std::vector<Kinesin*> > bound_tethered_;

		std::vector<event> events_;		// Holds all possible KMC events
		std::vector<int> kmc_list_;		// Holds which KMC functs to exe

	private:
		void GenerateMotors();
		void SetParameters();
		void InitializeLists();
		void InitializeEvents(); 

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
				system_properties *properties);

		int GetNumBoundUntethered();
		double GetWeight_BindTethered();
		double GetWeight_TetherBound();
		Kinesin* GetFreeMotor();
		Kinesin* GetBoundUntetheredMotor();

		void UpdateAllLists();
		void UpdateFreeTethered();
		void UpdateDocked();
		void UpdateDockedTethered();
		void UpdateBoundNULL();
		void UpdateBoundNULLTethered();
		void UpdateBoundATP();
		void UpdateBoundATP_Stalled();
		void UpdateBoundADPP_I();
		void UpdateBoundADPP_I_Stalled();
		void UpdateBoundADPP_I_Tethered();
		void UpdateBoundADPP_I_Tethered_Stalled();
		void UpdateBoundADPP_II();
		void UpdateBoundUntethered();
		void UpdateBoundTethered();

		void GenerateKMCList();
		void UpdateEvents();

		void RunKMC();
		void KMC_Bind_I();		// Bind free ADP head; convert to NULL
		void KMC_Bind_I_Tethered();
		void KMC_Bind_ATP();	// Bind ATP to NULL bound heads
		void KMC_Bind_ATP_Tethered(int x_dub); 
		void KMC_Hydrolyze(); 	// Convert ATP to ADPP on a bound head
		void KMC_Hydrolyze_Stalled();
		void KMC_Bind_II();	  // Bind docked head; other head must be ADPP
		void KMC_Bind_II_Tethered(int x_dub);
		void KMC_Unbind_II(); 		// Unbind ADPP heads; converts to ADP
		void KMC_Unbind_I();
		void KMC_Unbind_I_Stalled();
		void KMC_Unbind_I_Tethered(int x_dub);
		void KMC_Unbind_I_Tethered_Stalled(int x_dub);  //XXX
		void KMC_Tether_Free();
		void KMC_Tether_Bound();
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dub); 
};
#endif
