#include "binding_site.hpp"
#include "binding_head.hpp"
#include "catalytic_head.hpp"

bool BindingSite::HeadTrailing() { return occupant_->Trailing(); }