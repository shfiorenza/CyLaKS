#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"
#include <map>
#include <string>
#include <functional>

struct system_parameters;
struct system_properties;

// Make a structure to hold the number and type of population
struct pop_t{
	int n_entries_ = -1;
	std::string type_ = std::string("wut"); 
	int x_dist_ = -1;
	int x_dist_dub_ = -1;
}; 

class KinesinManagement{
	private:

	public:
		int n_motors_ = 0;
		
		// As a general rule, populations are 
		// untethered unless otherwise specified 
		int n_free_tethered_ = 0;
	 	int	n_bound_i_ = 0;
		int n_bound_i_bindable_ = 0; 
		int n_bound_ii_ = 0; 
		int n_bound_untethered_ = 0;
		int n_bound_ii_tethered_tot_ = 0;       //cut ?
		int n_stepable_ = 0; 
		// Tethered populations are organized by tether extension
		std::vector<int> n_bindable_to_teth_;   
		std::vector<int> n_bindable_fr_teth_; 
		std::vector<int> n_bound_i_tethered_;  	//cut ?
		std::vector<int> n_bound_ii_tethered_; 	//cut ?
		std::vector<int> n_bound_tethered_;	
		std::vector<int> n_stepable_to_teth_rest_;
		std::vector<int> n_stepable_fr_teth_rest_;

		// See kinesin header for description of below
		int dist_cutoff_;
		int comp_cutoff_;
		double rest_dist_;

		/* Diffusion stuff */
		double tau_;
		double p_diffuse_forward_; 
		double p_diffuse_backward_;
		// Diffusing away FROM rest means extension/compression increases,
		// whereas diffusing TO rest means extension/compression decreases
		std::vector<double> p_diffuse_to_teth_rest_;
		std::vector<double> p_diffuse_fr_teth_rest_;
	
		std::vector<int> diffusion_list_;

		/* Langmuir and kinematic stuff	*/
		double p_bind_i_;
		double p_bind_i_tethered_;
		double p_bind_ii_; 
		double p_unbind_i_; 
		double p_unbind_ii_;
		double p_tether_free_;
		double p_tether_bound_;
		double p_untether_free_;	
		double p_step_;
		// Indices refer to double the x_distance (in no. of sites)
		// between the stalk of a motor and the anchor of an xlink
		// e.g. p_untether_bound_[23] is for an x_dist of 11.5 sites
		std::vector<double> p_bind_ii_to_teth_;
		std::vector<double> p_bind_ii_fr_teth_;	
		std::vector<double> p_unbind_i_tethered_;  
		std::vector<double> p_unbind_ii_to_teth_; 
		std::vector<double> p_unbind_ii_fr_teth_;
		std::vector<double> p_untether_bound_; 	
		std::vector<double> p_step_to_teth_rest_;
		std::vector<double> p_step_fr_teth_rest_;

		std::vector<Kinesin> motors_; 
		std::vector<Kinesin*> free_tethered_list_;
		std::vector<Kinesin*> bound_i_list_; 
		std::vector<Kinesin*> bound_i_bindable_list_; 
		std::vector<Kinesin*> bound_ii_list_; 
		std::vector<Kinesin*> bound_untethered_;
		std::vector<Kinesin*> stepable_list_;
		// For our purposes, 'list' means one-dimensional vector, and
		// the 'table' means two-dimensions: tether extension and index
		std::vector<Kinesin*> bound_ii_tethered_list_; 
		std::vector< std::vector<Kinesin*> > bindable_to_teth_; 
		std::vector< std::vector<Kinesin*> > bindable_fr_teth_; 
		std::vector< std::vector<Kinesin*> > bound_i_tethered_;	// needed?
		std::vector< std::vector<Kinesin*> > bound_ii_tethered_table_;
		std::vector< std::vector<Kinesin*> > bound_tethered_; 
		std::vector< std::vector<Kinesin*> > stepable_to_rest_table_;
		std::vector< std::vector<Kinesin*> > stepable_fr_rest_table_;

		std::vector<int> kmc_list_;

		std::vector<pop_t> serial_pop_; 
		std::vector<pop_t> serial_kmc_;
		std::map<std::string, std::function<int(int)> > sampling_functs;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);
		void SetParameters();
		void GenerateMotors();
		void InitiateLists();
		void InitializeSerialPop(); 
		void InitializeFunctionMap();

		void UnboundCheck(Kinesin *motor);
		void BoundCheck(Kinesin *motor);

		int GetNumBoundUntethered();
		Kinesin* GetFreeMotor();
		Kinesin* GetBoundUntetheredMotor();

		void UpdateAllLists();
		void UpdateFreeTetheredList();
		void UpdateBoundIList(); 
		void UpdateBoundIBindableList();
		void UpdateBoundIIList();
		void UpdateBoundUntethered();
		void UpdateStepableList();
		void UpdateBindableToTeth(); 
		void UpdateBindableFromTeth(); 
		void UpdateBoundITethered(); 
		void UpdateBoundIITetheredList();
		void UpdateBoundIITetheredTable();
		void UpdateBoundTethered(); 
		void UpdateStepableTetheredTables();
		
		void GenerateDiffusionList();
		int GetNumToStepForward();
		int GetNumToStepBackward();
		int GetNumToStepToTethRest(int x_dist_doubled);
		int GetNumToStepFromTethRest(int x_dist_doubled); 

		void RunDiffusion();
		void RunDiffusion_Forward();
		void RunDiffusion_Backward();
		void RunDiffusion_To_Teth_Rest(int x_dist_doubled);
		void RunDiffusion_From_Teth_Rest(int x_dist_doubled);

		void GenerateKMCList();
		void UpdateSerializedPopulations(); 
		void UpdateSerializedEvents();
		double GetWeight_Bind_I_Tethered();
		double GetWeight_Tether_Bound();

		int GetNumToBind_I();
		int GetNumToBind_I_Tethered();
		int GetNumToBind_II();
		int GetNumToBind_II_To_Teth(int x_dist_doubled);   
		int GetNumToBind_II_From_Teth(int x_dist_doubled); 
		int GetNumToUnbind_I(); 
		int GetNumToUnbind_I_Tethered(int x_dist_doubled);
		int GetNumToUnbind_II();
		int GetNumToUnbind_II_To_Teth(int x_dist_doubled); 
		int GetNumToUnbind_II_From_Teth(int x_dist_doubled);
		int GetNumToTether_Free();
		int GetNumToTether_Bound(); 
		int GetNumToUntether_Free();
		int GetNumToUntether_Bound(int x_dist_doubled); 
		int GetNumToStep();
		int GetNumToStep_ToTethRest(int x_dist_doubled);
		int GetNumToStep_FromTethRest(int x_dist_doubled);

		void RunKMC();
		void KMC_Bind_I();
		void KMC_Bind_I_Tethered();
		void KMC_Bind_II();	
		void KMC_Bind_II_To_Teth(int x_dist_doubled);
		void KMC_Bind_II_From_Teth(int x_dist_doubled);	
		void KMC_Unbind_I();
		void KMC_Unbind_I_Tethered(int x_dist_doubled);
		void KMC_Unbind_II();
		void KMC_Unbind_II_To_Teth(int x_dist_doubled);	
		void KMC_Unbind_II_From_Teth(int x_dist_doubled);
		void KMC_Tether_Free();
		void KMC_Tether_Bound(); 
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dist_doubled); 
		void KMC_Step();
		void KMC_Step_ToTethRest(int x_dist_doubled);
		void KMC_Step_FromTethRest(int x_dist_doubled);
};
#endif
