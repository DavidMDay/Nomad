#pragma once

#include <vector>
struct TemporalParameters
{
  bool mConstraintErrorDiagnosticsFlag = false;
  bool mSkipSolve = false;
  bool mUseDefaultRho = true;
  double mRho = 0.0;
  double mStartTime = 0.0;
  int mConstraintCorrectionFrequency = 1;
  int mLineSkip = 0;
  int mPredictorCorrector = -1;
  static constexpr int default_nskips = 1;
  static constexpr int default_nsteps = 100;
  std::vector<double> mTimeStep;
  std::vector<int> mNumSkip;
  std::vector<int> mNumStep;
};
