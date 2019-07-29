#pragma once 
#include "event.h"
#include "entry.h"
struct system_parameters;
struct system_properties;
class Curator; 

class AssociatedProteinManagement{
	private:
		// Data types
		using POP_T = AssociatedProtein::Monomer;
		using ALT_T = Kinesin::head;
		using SITE_T = Tubulin;
		using MGMT_T = AssociatedProteinManagement;
		// ENTRY_T is defined in entry.h header
		using EVENT_T = Event<MGMT_T*, ENTRY_T>;
		// Use a template for 'Vec' rather than std::vector (aesthetic)
		template<class DATA_T> using Vec = std::vector<DATA_T>;
		// All possible KMC event objects; arbitrary sequential order
		Vec<EVENT_T> events_;
		// IDs of events; segregated by target pop. (for stat correction)
		Vec<Vec<int>> IDs_by_pop_;
		// Same as above but for secondary stat correction
		Vec<Vec<int>> IDs_by_root_;
		Vec<int*> n_avail_by_root_;
		// List of event IDs to execute any given timestep; dynamic
		Vec<int> IDs_to_exe_;
		// Temporarily holds entries after KMC events
		int n_scratched_ = 0;
		Vec<POP_T*> scratch_;
		
		// Pointers to global system params & props; same for all classes
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
		// WALLACE, MISTA
		Curator* wally = nullptr;

	public:
		int dist_cutoff_;		// see assoc. protein header
		int rest_dist_;			// see assoc. protein header
		int teth_cutoff_; 		// see kinesin header (dist_cutoff_ there)
		int comp_cutoff_;		// see kinesin header (comp_cutoff_ there)
		// Neighbor coop stuff; still kinda preliminary
		int max_neighbs_;
		double interaction_energy_;		// in kBT

		/* Population size trackers */ 
		int n_xlinks_ = 0; 		  // Total no. of xlink objects; static
		int n_active_ = 0; 		  // No. actively bound; dynamic
		int n_free_teth_ = 0;
		int n_bound_unteth_ = 0;
		int n_bound_i_teth_tot_ = 0;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		Vec<int> n_bound_i_;
		// Second index is [x] or [x_dub] for base or teth pops.
		Vec<Vec<int>> n_bound_i_teth_;
		Vec<Vec<int>> n_bound_ii_;
		// Second index is [x_dub], third & final index is [x]
		Vec<Vec<Vec<int>>> n_bound_ii_teth_same_; 
		Vec<Vec<Vec<int>>> n_bound_ii_teth_oppo_; 

		/* Probabilities of possible KMC events */
		double p_bind_i_teth_base_;
		double p_bind_ii_base_;
		double p_tether_free_; 
		double p_untether_free_;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		Vec<double> p_bind_i_;			
		Vec<double> p_unbind_i_;		
		Vec<double> p_diffuse_i_fwd_;
		Vec<double> p_diffuse_i_bck_;
		// Second index is [x] or [x_dub] for base or teth pops.
		Vec<Vec<double>> p_unbind_i_teth_;	
		Vec<Vec<double>> p_unbind_ii_;
		Vec<Vec<double>> p_diffuse_i_to_teth_rest_;
		Vec<Vec<double>> p_diffuse_i_fr_teth_rest_;
		Vec<Vec<double>> p_diffuse_ii_to_rest_;
		Vec<Vec<double>> p_diffuse_ii_fr_rest_;
		// Second index is [x_dub]; third & final index is [x]
		Vec<Vec<Vec<double>>> p_unbind_ii_to_teth_;
		Vec<Vec<Vec<double>>> p_unbind_ii_fr_teth_;
		Vec<Vec<Vec<double>>> p_diffuse_ii_to_both_;
		Vec<Vec<Vec<double>>> p_diffuse_ii_fr_both_;
		Vec<Vec<Vec<double>>> p_diffuse_ii_to_self_fr_teth_;
		Vec<Vec<Vec<double>>> p_diffuse_ii_fr_self_to_teth_;

		/* Lists that track different population types */
		Vec<AssociatedProtein> xlinks_;				// Actual xlink objects
		Vec<AssociatedProtein*> active_;
		Vec<AssociatedProtein*> bound_unteth_;
		Vec<ENTRY_T> free_teth_;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		// Second index is actual xlink entry
		Vec<Vec<ENTRY_T>> bound_i_;
		// Second index is [x] or [x_dub]; third index is xlink entry
		Vec<Vec<Vec<ENTRY_T>>> bound_ii_;
		Vec<Vec<Vec<ENTRY_T>>> bound_i_teth_;
		// Second index is [x_dub]; third index is [x]; fourth is entry
		// e.g., [0][16][2][1] -> 2nd xlink w/ x_dub=16, x=2, & 0 neighbs
		Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_oppo_;
		Vec<Vec<Vec<Vec<ENTRY_T>>>> bound_ii_teth_same_;

	private:
		void GenerateXLinks();
		void SetParameters();
		void InitializeLists();
		void InitializeEvents();

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		AssociatedProtein* GetFreeXlink();
		AssociatedProtein* GetBoundUntetheredXlink();
		double GetWeight_Bind_II(); 
		double GetWeight_Bind_I_Teth();
		double GetWeight_Bind_II_Teth();

		POP_T* CheckScratchFor(std::string pop);
		void SaveToScratch(POP_T* head);

		void Update_All_Lists();
		void Update_List_Relay(std::string event, std::string target_pop);
		void Update_Free_Teth();
		void Update_Bound_Unteth(); 
		void Update_Bound_I();
		void Update_Bound_I_Teth();
		void Update_Bound_II();
		void Update_Bound_II_Teth();
		void Update_Bind_II_Candidate();
		void Update_Bind_I_Teth_Candidate();
		void Update_Bind_II_Teth_Candidate();	//XXX add

		void Run_KMC();
		void Refresh_Populations(); 
		void Generate_Execution_Sequence();
		int Sample_Event_Statistics(); 

		void Execute_Function_Relay(ENTRY_T target, int code); 
		void Diffuse(POP_T* head, int dir);
		void Bind_I(SITE_T* target_site);
		void Bind_II(POP_T* bound_head);
		void Unbind_I(POP_T* bound_head);
		void Unbind_II(POP_T* bound_head);
		void Bind_I_Teth(POP_T* satellite_head);
		void Tether_Free(ALT_T* untethered_head); 
		void Untether_Free(POP_T* satellite_head);
};
