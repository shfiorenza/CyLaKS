#include "cylaks/filament_tester.hpp"
#include "cylaks/protein_tester.hpp"

void FilamentTester::Initialize(ProteinTester *proteins) {

  proteins_ = proteins;
  FilamentManager::proteins_ = dynamic_cast<ProteinManager *>(proteins_);
  FilamentManager::SetParameters();
  FilamentManager::GenerateFilaments();
}
