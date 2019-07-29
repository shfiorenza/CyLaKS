#pragma once
#include "tubulin.h"
#include <variant>

using 
ENTRY_T = std::variant<Tubulin*, AssociatedProtein::Monomer*, Kinesin::head*>;
