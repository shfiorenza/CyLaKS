#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"
#include <map>
#include <string>
#include <functional>

struct system_parameters;
struct system_properties;

//XXX to-do:
//XXX   - make UpdateAllLists() one for loop w/ a shitton of conditions
//XXX	- experiment with #pragma omp parallel for in UpdateLists()
//XXX   - investigate using sampling_functs_['name'] rather than find()
//XXX 	- reorganize KMC functions so non-tether group is first

class KinesinManagement{
	private:
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
		// Structure to hold the size and type of population
		struct pop_t{
			int n_entries_ = -1;
			std::string type_ = std::string("wut"); 
			int x_dist_ = -1;
			int x_dist_dub_ = -1;
		}; 

	public:
		int n_motors_ = 0;
		// motors actively bound to some MT or xlink; dynamically updated
		int n_active_ = 0;	
		
		// As a general rule, populations are 
		// untethered unless otherwise specified 
		int n_free_tethered_ = 0;
	 	int	n_bound_i_ = 0;
		int n_bound_i_bindable_ = 0; 
		int n_bound_ii_ = 0; 
		int n_bound_untethered_ = 0;
		int n_stepable_ = 0; 
		// Indices refer to double the number of sites (x_dist_doubled_)
		// between the stalk of a motor and the anchor of an xlink
		// e.g. n_bound_tethered_[23] is for a distance of 11.5 sites
		std::vector<int> n_bound_i_tethered_;
		std::vector<int> n_bound_i_bindable_to_teth_;   
		std::vector<int> n_bound_i_bindable_fr_teth_; 
		std::vector<int> n_bound_ii_tethered_;
		std::vector<int> n_bound_tethered_;
		std::vector<int> n_stepable_to_teth_;
		std::vector<int> n_stepable_fr_teth_;

		// See kinesin header for description of below
		int dist_cutoff_;
		int comp_cutoff_;
		double rest_dist_;

		// Kinematics probabilities 
		double p_bind_i_;
		double p_bind_i_tethered_;
		double p_bind_ii_; 
		double p_unbind_i_; 
		double p_unbind_ii_;
		double p_tether_free_;
		double p_tether_bound_;
		double p_untether_free_;	
		double p_step_;
		// KMC events involving tethers have diff. rates based on extension
		std::vector<double> p_bind_ii_to_teth_;
		std::vector<double> p_bind_ii_fr_teth_;	
		std::vector<double> p_unbind_i_tethered_;  
		std::vector<double> p_unbind_ii_to_teth_; 
		std::vector<double> p_unbind_ii_fr_teth_;
		std::vector<double> p_untether_bound_; 	
		std::vector<double> p_step_to_teth_;
		std::vector<double> p_step_fr_teth_;

		// 1-D vectors, index is simply motor entry
		std::vector<Kinesin> motors_; 
		std::vector<Kinesin*> active_;
		std::vector<Kinesin*> free_tethered_;
		std::vector<Kinesin*> bound_i_; 
		std::vector<Kinesin*> bound_i_bindable_; 
		std::vector<Kinesin*> bound_ii_; 
		std::vector<Kinesin*> bound_untethered_;
		std::vector<Kinesin*> stepable_;
		// 2-D vectors: first index is x_dist_doubled, second is motor entry
		// e.g. bound_i_tethered_[16][0] is the 1st motor with 2x = 16
		std::vector< std::vector<Kinesin*> > bound_i_tethered_;
		std::vector< std::vector<Kinesin*> > bound_i_bindable_to_teth_; 
		std::vector< std::vector<Kinesin*> > bound_i_bindable_fr_teth_; 
		std::vector< std::vector<Kinesin*> > bound_ii_tethered_;
		std::vector< std::vector<Kinesin*> > bound_tethered_; 
		std::vector< std::vector<Kinesin*> > stepable_to_teth_;
		std::vector< std::vector<Kinesin*> > stepable_fr_teth_;

		// Holds KMC events encoded in integers
		std::vector<int> kmc_list_;
		// Serialized (in regards to teth extension) vectors that store
		// number of entries, the population/event label, and x_dub
		// (facilitates the parallelization of statistical sampling) 
		std::vector<pop_t> serial_pop_; 
		std::vector<pop_t> serial_kmc_;
		// Map of population/event label to sampling function;
		// see InitializeFUnctionMap() for function definitions
		std::map<std::string, std::function<int(int)> > sampling_functs_;

	private:
		void GenerateMotors();
		void SetParameters();
		void InitializeLists();
		void InitializeSerialPop(); 
		void InitializeSamplingFunctions();

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void UnboundCheck(Kinesin *motor);
		void BoundCheck(Kinesin *motor);

		int GetNumBoundUntethered();
		Kinesin* GetFreeMotor();
		Kinesin* GetBoundUntetheredMotor();

		void UpdateAllLists();
		void UpdateFreeTethered();
		void UpdateBoundI(); 
		void UpdateBoundIBindable();
		void UpdateBoundII();
		void UpdateBoundUntethered();
		void UpdateStepable();
		void UpdateBindableToTeth(); 
		void UpdateBindableFromTeth(); 
		void UpdateBoundITethered(); 
		void UpdateBoundIITethered();
		void UpdateBoundTethered(); 
		void UpdateStepableTethered();

		void GenerateKMCList();
		void UpdateSerializedPopulations(); 
		void UpdateSerializedEvents();
		double GetWeightBindITethered();
		double GetWeightTetherBound();

		void RunKMC();
		void KMC_Bind_I();
		void KMC_Bind_I_Tethered();
		void KMC_Bind_II();	
		void KMC_Bind_II_To_Teth_Rest(int x_dist_doubled);
		void KMC_Bind_II_From_Teth_Rest(int x_dist_doubled);	
		void KMC_Unbind_I();
		void KMC_Unbind_I_Tethered(int x_dist_doubled);
		void KMC_Unbind_II();
		void KMC_Unbind_II_To_Teth_Rest(int x_dist_doubled);	
		void KMC_Unbind_II_From_Teth_Rest(int x_dist_doubled);
		void KMC_Tether_Free();
		void KMC_Tether_Bound(); 
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dist_doubled); 
		void KMC_Step();
		void KMC_Step_To_Teth_Rest(int x_dist_doubled);
		void KMC_Step_From_Teth_Rest(int x_dist_doubled);
};
#endif
