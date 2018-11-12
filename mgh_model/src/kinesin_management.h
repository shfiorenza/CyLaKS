#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"

struct system_parameters;
struct system_properties;

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
		int n_stalled_ = 0;						//cut 
		int n_stepable_ = 0; 
		// Tethered populations are organized by tether extension
		std::vector<int> n_bindable_to_teth_;   
		std::vector<int> n_bindable_from_teth_; 
		std::vector<int> n_bound_i_tethered_;  	//cut ?
		std::vector<int> n_bound_ii_tethered_; 	//cut ?
		std::vector<int> n_bound_tethered_;	
		std::vector<int> n_stepable_to_teth_rest_;
		std::vector<int> n_stepable_from_teth_rest_;
		std::vector<int> n_stalled_to_teth_rest_;	// cut
		std::vector<int> n_stalled_from_teth_rest_;	// cut

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
		std::vector<double> p_diffuse_from_teth_rest_;
	
		std::vector<int> diffusion_list_;

		/* Langmuir and kinematic stuff	*/
		double p_bind_i_;
		double p_bind_i_tethered_;
		double p_bind_ii_; 
		std::vector<double> p_bind_ii_to_teth_;
		std::vector<double> p_bind_ii_from_teth_;	
		double p_unbind_i_; 
		std::vector<double> p_unbind_i_tethered_;  
		double p_unbind_ii_stepable_ ; 				 // collapse into 1 <
		double p_unbind_ii_stalled_; 			
		std::vector<double> p_unbind_ii_to_teth_; 
		std::vector<double> p_unbind_ii_from_teth_;
		double p_tether_free_;
		double p_tether_bound_;
		// Indices refer to double the x_distance (in no. of sites)
		// between the stalk of a motor and the anchor of an xlink
		// e.g. p_untether_bound_[23] is for an x_dist of 11.5 sites
		double p_untether_free_;	
		std::vector<double> p_untether_bound_; 	
		double p_step_;
		double p_failstep_;							// cut
		std::vector<double> p_step_to_teth_rest_;
		std::vector<double> p_step_from_teth_rest_;
		std::vector<double> p_failstep_to_teth_rest_; // cut
		std::vector<double> p_failstep_from_teth_rest_; // cut

		std::vector<Kinesin> motors_; 
		std::vector<Kinesin*> free_tethered_list_;
		std::vector<Kinesin*> bound_i_list_; 
		std::vector<Kinesin*> bound_i_bindable_list_; 
		std::vector<Kinesin*> bound_ii_list_; 
		std::vector<Kinesin*> bound_untethered_;
		std::vector<Kinesin*> stepable_list_;
		std::vector<Kinesin*> stalled_list_;		// cut
		// For our purposes, 'list' means one-dimensional vector, and
		// the 'table' means two-dimensions: tether extension and index
		std::vector<Kinesin*> bound_ii_tethered_list_; 
		std::vector< std::vector<Kinesin*> > bindable_to_teth_; 
		std::vector< std::vector<Kinesin*> > bindable_from_teth_; 
		std::vector< std::vector<Kinesin*> > bound_i_tethered_;	// needed?
		std::vector< std::vector<Kinesin*> > bound_ii_tethered_table_;
		std::vector< std::vector<Kinesin*> > bound_tethered_; 
		std::vector< std::vector<Kinesin*> > stepable_to_rest_table_;
		std::vector< std::vector<Kinesin*> > stepable_from_rest_table_;
		std::vector< std::vector<Kinesin*> > stalled_to_rest_table_; // cut
		std::vector< std::vector<Kinesin*> > stalled_from_rest_table_; // cut

		std::vector<int> kmc_list_;

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

		void UnboundCheck(Kinesin *motor);
		void BoundCheck(Kinesin *motor);

		int GetNumBoundUntethered();
		Kinesin* GetFreeMotor();
		Kinesin* GetBoundUntetheredMotor();

		void UpdateFreeTetheredList();
		void UpdateBoundIList(); 
		void UpdateBoundIBindableList();
		void UpdateBoundIIList();
		void UpdateBoundUntethered();
		void UpdateStepableList();
		void UpdateStalledList();	 
		void UpdateBindableToTeth(); 
		void UpdateBindableFromTeth(); 
		void UpdateBoundITethered(); 
		void UpdateBoundIITetheredList();
		void UpdateBoundIITetheredTable();
		void UpdateBoundTethered(); 
		void UpdateStepableTetheredTables();
		void UpdateStalledTetheredTables();
		
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
		// Roman numerals refer to binding stage being executed 
		int GetNumToBind_I();
		int GetNumToBind_I_Tethered();
		int GetNumToBind_II();
		int GetNumToBind_II_To_Teth(int x_dist_doubled);   
		int GetNumToBind_II_From_Teth(int x_dist_doubled); 
		int GetNumToUnbind_I(); 
		int GetNumToUnbind_I_Tethered(int x_dist_doubled);
		int GetNumToUnbind_II_Stepable();
		int GetNumToUnbind_II_Stalled();
		int GetNumToUnbind_II_To_Teth(int x_dist_doubled); 
		int GetNumToUnbind_II_From_Teth(int x_dist_doubled);
		int GetNumToTether_Free();
		int GetNumToTether_Bound(); //XXX
		int GetNumToUntether_Free();
		int GetNumToUntether_Bound(int x_dist_doubled); //XXX
		int GetNumToStep();
		int GetNumToFailstep();  //cut
		int GetNumToStep_ToTethRest(int x_dist_doubled);
		int GetNumToStep_FromTethRest(int x_dist_doubled);
		int GetNumToFailstep_ToTethRest(int x_dist_doubled); //cut
		int GetNumToFailstep_FromTethRest(int x_dist_doubled); //cut

		void RunKMC();
		void KMC_Bind_I();
		void KMC_Bind_I_Tethered();
		void KMC_Bind_II();		// FIXME to only deal w/ untethered binds
		void KMC_Bind_II_To_Teth(int x_dist_doubled);		//XXX
		void KMC_Bind_II_From_Teth(int x_dist_doubled);		//XXX
		void KMC_Unbind_I();	// FIXME to only deal w/ untethered unbinds
		void KMC_Unbind_I_Tethered(int x_dist_doubled);		//XXX
		void KMC_Unbind_II_Stepable();
		void KMC_Unbind_II_Stalled(); 
		void KMC_Unbind_II_To_Teth(int x_dist_doubled);		//XXX
		void KMC_Unbind_II_From_Teth(int x_dist_doubled);	//XXX
		void KMC_Tether_Free();
		void KMC_Tether_Bound(); 					 //XXX
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dist_doubled); //XXX
		void KMC_Step();
		void KMC_Failstep();								//cut
		void KMC_Step_ToTethRest(int x_dist_doubled);
		void KMC_Step_FromTethRest(int x_dist_doubled);
		void KMC_Failstep_ToTethRest(int x_dist_doubled);	//cut
		void KMC_Failstep_FromTethRest(int x_dist_doubled); //cut
};
#endif
