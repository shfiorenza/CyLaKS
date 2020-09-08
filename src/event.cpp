#include "event.hpp"
#include "associated_protein.h"
#include "kinesin.h"
#include "tubulin.h"

template class Event<
    std::variant<Tubulin *, AssociatedProtein::Monomer *, Kinesin::Monomer *>>;
