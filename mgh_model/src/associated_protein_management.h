#pragma once 
#include <string>
#include <variant>
#include <functional>
#include "event.h"
#include "associated_protein.h"
struct system_parameters;
struct system_properties;
class Curator; 

class AssociatedProteinManagement{
	private:
		// Aliases to make our work a little more neat
		using POP_T = AssociatedProtein::Monomer;
		using SITE_T = Tubulin;
		using MGMT_T = AssociatedProteinManagement;
		using EVENT_T = Event<POP_T*, SITE_T*, MGMT_T*>;
		using Entry = EVENT_T::Entry;
		template<class DATA_T> using Vec = std::vector<DATA_T>;
		// All possible KMC event objects; arbitrary sequential order
		Vec<EVENT_T> events_;
		// IDs of events; segregated by target population
		Vec<Vec<int>> IDs_by_pop_;
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
		int n_heads_i_teth_tot_ = 0;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		Vec<int> n_heads_i_;
		// Second index is [x] or [x_dub] for base or teth pops.
		Vec<Vec<int>> n_heads_i_teth_;
		Vec<Vec<int>> n_heads_ii_;
		// Second index is [x_dub], third & final index is [x]
		Vec<Vec<Vec<int>>> n_heads_ii_teth_same_; 
		Vec<Vec<Vec<int>>> n_heads_ii_teth_oppo_; 

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
		Vec<AssociatedProtein*> free_teth_;
		Vec<AssociatedProtein*> bound_unteth_;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		// Second index is actual xlink entry
		Vec<Vec<Entry>> heads_i_;
		// Second index is [x] or [x_dub]; third index is xlink entry
		Vec<Vec<Vec<POP_T*>>> heads_ii_;
		Vec<Vec<Vec<POP_T*>>> heads_i_teth_;
		// Second index is [x_dub]; third index is [x]; fourth is entry
		// e.g., [0][16][2][1] -> 2nd xlink w/ x_dub=16, x=2, & 0 neighbs
		Vec<Vec<Vec<Vec<POP_T*>>>> heads_ii_teth_oppo_;
		Vec<Vec<Vec<Vec<POP_T*>>>> heads_ii_teth_same_;

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

		void Update_Relay(std::string pop);
		void Update_All();
		void Update_Free_Teth();
		void Update_Bound_Unteth(); 
		void Update_Heads_I();
		void Update_Heads_I_Teth();
		void Update_Heads_II();
		void Update_Heads_II_Teth();

		void Run_KMC();
		void Refresh_Population_Labels(); 
		void Generate_Execution_Sequence();
		void KMC_Relay(Entry target, int code); 
		void SaveToScratch(POP_T* head);
		// Diffusion events
		void Diffuse_I(POP_T* head, int dir);
		void Diffuse_II_To_Rest(int n_neighbs, int x);
		void Diffuse_II_Fr_Rest(int n_neighbs, int x);
		void Diffuse_I_To_Teth(int n_neighbs, int x_dub);
		void Diffuse_I_Fr_Teth(int n_neighbs, int x_dub);
		void Diffuse_II_To_Both(int n_neighbs, int x_dub, int x);
		void Diffuse_II_Fr_Both(int n_neighbs, int x_dub, int x);
		void Diffuse_II_To_Self_Fr_Teth(int n_neighbs, int x_dub, int x);
		void Diffuse_II_Fr_Self_To_Teth(int n_neighbs, int x_dub, int x);
		// Kinematic events
		void Bind_I(SITE_T* site);
		void Bind_II();		//XXX add neighb coop in weights
		void Unbind_I(POP_T* head);
		void Unbind_II(int n_neighbs, int x);
		void Bind_I_Teth(); //XXX add neighb coop in weights
		void Bind_II_Teth(); //XXX add neighb coop in weights
		void Unbind_I_Teth(int n_neighbs, int x_dub); 
		void Unbind_II_To_Teth(int n_neighbs, int x_dub, int x);
		void Unbind_II_Fr_Teth(int n_neighbs, int x_dub, int x);
		void Tether_Free(); 
		void Untether_Free();
};
