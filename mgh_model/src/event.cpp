#include "event.hpp"
#include "associated_protein_management.h"
#include "kinesin_management.h"

template class Event<
    AssociatedProteinManagement *,
    std::variant<Tubulin *, AssociatedProtein::Monomer *, Kinesin::head *>>;

template class Event<
    KinesinManagement *,
    std::variant<Tubulin *, AssociatedProtein::Monomer *, Kinesin::head *>>;