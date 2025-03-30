#include <array>
#include <vector>
#include <iostream>
#include "FourierCoefficients.hpp"
#include "TemporalParameters.hpp"
#include "FixedStepIntegrator.hpp"
#include "Tests.hpp"
#include "SpaceTimeTransient.hpp"
#include "HB.hpp"
#include <cassert>

TemporalParameters get_one_fourth_slow(){
  TemporalParameters t;
  constexpr int n = 128;
  // omega = 1/4,  t/4 = 2 pi
  constexpr int num_cycles = 100;
  const double h = 8.0 * M_PI/static_cast<double>(n);
  t.mTimeStep.emplace_back(h);
  t.mNumSkip.emplace_back(1);
  t.mNumStep.emplace_back(n*num_cycles);
  return t;
}

TemporalParameters get_one_fourth_tp(){
  TemporalParameters t;
  constexpr int n = 128;
  // omega = 1/4,  t/4 = 2 pi
  constexpr int num_cycles = 2;
  const double h = 8.0 * M_PI/static_cast<double>(n);
  t.mTimeStep.emplace_back(h);
  t.mNumSkip.emplace_back(1);
  t.mNumStep.emplace_back(n*num_cycles);
  return t;
}

TemporalParameters get_one_third_slow(){
  TemporalParameters t;
  constexpr int n = 128;
  // omega = 1/3,  t/3 = 2 pi
  constexpr int num_cycles = 100;
  const double h = 6.0 * M_PI/static_cast<double>(n);
  t.mTimeStep.emplace_back(h);
  t.mNumSkip.emplace_back(1);
  t.mNumStep.emplace_back(n*num_cycles);
  return t;
}

TemporalParameters get_one_third_tp(){
  TemporalParameters t;
  constexpr int n = 128;
  // omega = 1/3,  t/3 = 2 pi
  constexpr int num_cycles = 2;
  const double h = 6.0 * M_PI/static_cast<double>(n);
  t.mTimeStep.emplace_back(h);
  t.mNumSkip.emplace_back(1);
  t.mNumStep.emplace_back(n*num_cycles);
  return t;
}

int TestLiuOde() {
  auto param = get_one_fourth_tp();// omega = 1/4
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing(); //u'' + 2 xi u' + u + u^3 = f
  xx.circular_frequency = 0.25;
  // 100 cycles,   to = 6.49378,  uo = 0.0266048
  std::vector<double> u0 = {0.0};
  std::vector<double> v0 = {0.0};
  Ode.SetInitialCondition(u0, v0);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  return EXIT_SUCCESS;
}

int TestLightlyDampedOde() {
  auto param = get_one_fourth_tp();// omega = 1/4
  param.mNumStep[0] = 128; // 1 cycle
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing_low(); //u'' + 2 xi u' + u + u^3 = f
  xx.circular_frequency = 0.25;
  param.mStartTime = 6.30261;
  std::vector<double> u0 = {0.026641};
  std::vector<double> v0 = {0.0};
  Ode.SetInitialCondition(u0, v0);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  return EXIT_SUCCESS;
}

int TestSpaceTimeTransOde() {
  auto param = get_one_fourth_tp();// omega = 1/4
  param.mNumStep[0] = 128; // 1 cycle
  param.mStartTime = 6.49378;
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing(); //u'' + 2 xi u' + u + u^3 = f
  xx.circular_frequency = 0.25;
  std::vector<double> u0 = {0.0266048};
  std::vector<double> v0 = {0.0};
  Ode.SetInitialCondition(u0, v0);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  return EXIT_SUCCESS;
}

// xi = 1/100
// w=1/4, 8 pi,  6.3,  6.3/8 ~ 0.8  to = 6.30261; u0 = {0.026641};
// w=1/3, 6 pi,  4.6   4.6/6 ~ 0.8
// xi = 1/10                        to = 6.49378; uo = {0.0266048};

int TestIncreasedOmega() {
  auto param = get_one_third_tp();// omega = 1/3
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing_low(); //u'' + 2 xi u' + u + u^3 = f
  xx.circular_frequency = 1.0/3.0;
  // 100 cycles,   to = 4.60127, uo = 0.0281483
  std::vector<double> u0 = {0.0};
  std::vector<double> v0 = {0.0};
  Ode.SetInitialCondition(u0, v0);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  return EXIT_SUCCESS;
}


int TestIncreasedOmegaLiu() {
  auto param = get_one_third_tp();// omega = 1/3
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing(); //u'' + 2 xi u' + u + u^3 = f
  xx.circular_frequency = 1.0/3.0; // (to,uo) = 4.91505 , 0.0280293
  std::vector<double> u0 = {0.0};
  std::vector<double> v0 = {0.0};
  Ode.SetInitialCondition(u0, v0);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  return EXIT_SUCCESS;
}
// I learned that Ode::Force(step) returns the
// force at step+1/2.  ... There needs to be
// two different objects,  Force(step) for Harmonic
// Balance and force(step+1/2) for Trapezoid 
// ... the time marching algorithm for
// finding solution trajectories.
int TestLinearized() {
  auto param = get_one_fourth_tp();// omega = 1/4
  int neq = 1;
  SpaceTimeTransient Ode = SpaceTimeTransient(param, neq);
  HB xx = get_duffing_linear(); //u'' + 2 xi u' + u = f
  assert( xx.linear );
  duffing_linear(xx);
  std::vector<double> u0 = {-0.00141819};
  std::vector<double> v0 = {0.00664776};
  Ode.SetInitialCondition(u0, v0);
  //Ode.OutputTrajectory(true);
  Ode.Initialize(xx);
  Ode.StoreSolution();
  Ode.Solve(xx);
  for (int order = 1; order < 5; order++)
  {
    std::cout << "order " << order << "   ";
    xx.order = order;
    Ode.Wip(xx);
  }
  return EXIT_SUCCESS;
}

