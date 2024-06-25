#ifndef _CYLAKS_PROTEIN_TESTER_
#define _CYLAKS_PROTEIN_TESTER_
#include "protein_manager.hpp"

class FilamentTester;

class ProteinTester : public ProteinManager {
protected:
  Map<Str, Vec<Pair<size_t, size_t>>> test_stats_;
  Map<Str, Vec<double>> test_ref_;

  FilamentTester *filaments_{nullptr};

public:
  // Reservoir<TestMotor> motors_;

protected:
  void UpdateFilaments();
  void ReportTestStatistics();

  void SetTestMode();
  void InitializeTest_Filament_Ablation();
  void InitializeTest_Filament_Separation();
  void InitializeTest_Filament_ForcedSlide();
  void InitializeTest_Filament_HeteroTubulin();
  void InitializeTest_Motor_Heterodimer();
  void InitializeTest_Motor_LatticeStep();
  void InitializeTest_Motor_LatticeBind();
  void InitializeTest_Xlink_Diffusion();
  void InitializeTest_Xlink_Diffusion_Boltzmann();
  void InitializeTest_Xlink_Bind_II();
  void InitializeTest_Shepherding();

public:
  ProteinTester() : ProteinManager() {}
  ~ProteinTester() { ReportTestStatistics(); }
  void Initialize(FilamentTester *filaments);
};
#endif