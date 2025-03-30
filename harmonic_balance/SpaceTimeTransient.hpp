#pragma once

#include <valarray>
#include "TemporalParameters.hpp"
#include "FixedStepIntegrator.hpp"
#include "HB.hpp"
#include <vector>

class SpaceTimeTransient {
 public:
  explicit SpaceTimeTransient(TemporalParameters p, int neq = 1);
  virtual ~SpaceTimeTransient() = default;

  double Force(int step, HB x); // Force step+1/2


  void PrintMe();
  std::vector<double> Load(int step, HB x);

// ODE + Trapezoid rule tasks
  std::vector<double> OdeMatrix(double xi, double u);
  std::vector<double> TrapezoidRhs(HB xx); // rhs = yo + (h/2)(A(yo)yo + f0 + f1)
  std::vector<double> NonlinearResidual(HB xx);  //(I - (h/2) A(y1)) y1 - rhs
  std::vector<double> JacobianMatrix(HB xx);// analytic
  std::vector<double> fdJacobianMatrix(HB xx, double tau);// finite difference
  int TrapezoidPredictor();
  int Solve(HB xx);  // Solve Ode
  int NonlinearSolveAtEachStep(HB xx);
  void find_max_disp( int step);

  int Wip(HB xx);

  void SetInitialCondition( std::vector<double>& deformation, std::vector<double>& velocity);
  int Initialize(HB xx);
 
  void StoreSolution() ;
  void OutputTrajectory(bool default_is_false);
  int write_trajectory_to_disk(double step_size);
  void Print4Vectors( const std::vector<double>& t,
    const std::vector<double>& d,
    const std::vector<double>& v, const std::vector<double>& a);
  bool write_to_file = false;
  bool hasdamping = false;
  FixedStepIntegrator mTime_integ;
  int numIntegrators = 0;
  int mNeq = 1;
  int mNumCycles = 1;
  std::valarray<double> mSpaceTimeDisp;
  std::valarray<double> mSpaceTimeVel;
  std::valarray<double> mSpaceTimeAccel;
  std::valarray<double> mSpaceTimeForce;
  std::vector<double> mD_n;
  std::vector<double> mVel_n;
  std::vector<double> mAccel_n;
  std::vector<double> mAccel_n_post;
  std::vector<double> mASetForce;
};



