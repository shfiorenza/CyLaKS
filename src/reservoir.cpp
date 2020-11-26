#include "reservoir.hpp"

template class Reservoir<Protein>;
template class Reservoir<Motor>;

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SetParameters() {}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SetWeights() {}