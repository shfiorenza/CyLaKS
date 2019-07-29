#include "event.h"
#include "associated_protein_management.h" 

template class
Event<AssociatedProteinManagement*, 
	  std::variant<Tubulin*, AssociatedProtein::Monomer*, Kinesin::head*>>;
