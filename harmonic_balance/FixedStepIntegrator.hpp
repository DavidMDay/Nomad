#pragma once

#include <array>   // for array
#include <vector>  // for vector

#include "TemporalParameters.hpp"

struct FixedStepIntegrator {
 public:
  FixedStepIntegrator() = default;
  FixedStepIntegrator(TemporalParameters t);
  virtual ~FixedStepIntegrator() = default;
  double InitialStepSize() const { return mInitialStepSize; }
  double StepSize() const { return mCurrentStepSize; }
  double TotalDuration() const { return static_cast<double>(mNsteps) * mInitialStepSize; }
  double TimeAtStartOfStep() const { return mTime; }
  int NSteps() const { return mNsteps; }
  int StepsTaken() const { return mStepsTaken; }
  int TotalStepsTaken() const { return mTotalStepsTaken; }
  std::vector<double> getTimeVector() const;
  int timeVectorSize() const;
  void Reset();
  void Set(TemporalParameters t);
  void IncrementTimeBeforeCallingTakeNextStep();
  void TakeNextStep();
  void PrepareLastStep(double timeAtEndOfStep, double initialTime) const;
  void PrintMe() const;

  double mTime = 0;
  double mCurrentStepSize = 0.0;
  double mInitialStepSize = 0.0;
  int mNsteps = 0;
  int mStepsTaken = 0;
  int mTotalStepsTaken = 0;
  int mNumSkip = 0;
};
