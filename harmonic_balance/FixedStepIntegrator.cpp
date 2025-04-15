#include "FixedStepIntegrator.hpp"
#include "TemporalParameters.hpp"
#include <iostream>

FixedStepIntegrator::FixedStepIntegrator(TemporalParameters t)
    : mTime(t.mStartTime),
      mCurrentStepSize(t.mTimeStep[0]),
      mInitialStepSize(t.mTimeStep[0]),
      mNsteps(t.mNumStep[0]),
      mNumSkip(t.mNumSkip[0]) {}

void FixedStepIntegrator::Set(TemporalParameters t) {
  mTime = t.mStartTime;
  mCurrentStepSize = t.mTimeStep[0];
  mInitialStepSize = t.mTimeStep[0];
  mNsteps = t.mNumStep[0];
  mNumSkip = t.mNumSkip[0];
}

void FixedStepIntegrator::Reset() {
  mCurrentStepSize = mInitialStepSize;
  mStepsTaken = 0;
  mTotalStepsTaken = 0;
}

std::vector<double> FixedStepIntegrator::getTimeVector() const {
  int period = (mNumSkip >= mNsteps ? 1 : mNumSkip);
  int m = mNsteps / period;
  std::vector<double> time(m, 0);
  for (int k = 0; k < m; k++) {
    time[k] = static_cast<double>(period * (k + 1)) * mInitialStepSize;
  }
  return time;
}

int FixedStepIntegrator::timeVectorSize() const {
  int period = (mNumSkip >= mNsteps ? 1 : mNumSkip);
  return mNsteps / period;
}

void FixedStepIntegrator::IncrementTimeBeforeCallingTakeNextStep() { mTime += StepSize(); }

void FixedStepIntegrator::TakeNextStep() {
  mStepsTaken++;
  mTotalStepsTaken++;
  if (mNsteps == mStepsTaken) {
    mCurrentStepSize = 0.0;
    std::cout << " FixedStepIntegrator: this one is fine but can no longer take a step\n";
  }
}

void FixedStepIntegrator::PrepareLastStep(double timeAtEndOfStep, double initialTime) const {
  double targetTime = initialTime + TotalDuration();
  if (timeAtEndOfStep >= targetTime) {
    std::cout << " FixedStepIntegrator:set last step missing\n";
  }
}

void FixedStepIntegrator::PrintMe() const {
  std::cout << "FixedStepIntegrator: t = " << mTime << " step " << mStepsTaken << " of " << mNsteps
            << "\n";
  std::cout << "      step sizes " << mCurrentStepSize << ", " << mInitialStepSize << " total "
            << mTotalStepsTaken << " skip " << mNumSkip << "\n";
}