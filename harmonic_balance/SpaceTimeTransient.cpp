#include "SpaceTimeTransient.hpp"

#include <cstdio>
#include <cstdlib>

#include <QRdecomposition.hpp>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <triuSol.hpp>
#include <util.hpp>
#include <utility>

#include "HB.hpp"
#include "Operators.hpp"
#include "matrix2x2.hpp"

SpaceTimeTransient::SpaceTimeTransient(TemporalParameters param, int neq) : mNeq(neq) {
  mTime_integ.Set(param);
  try {
    mSpaceTimeDisp.resize(mTime_integ.NSteps() * mNeq, 0.0);
    mSpaceTimeVel.resize(mTime_integ.NSteps() * mNeq, 0.0);
    mSpaceTimeAccel.resize(mTime_integ.NSteps() * mNeq, 0.0);
    mSpaceTimeForce.resize(mTime_integ.NSteps() * mNeq, 0.0);
  } catch (std::bad_alloc &own) {
    std::cout << "SpaceTimeTransient: insufficient memory " << own.what() << "\n";
  } catch (std::exception &hah) {
    std::cout << "SpaceTimeTransient: caught exception " << hah.what() << "\n";
  }
}

int SpaceTimeTransient::Initialize(HB xx) {
  mASetForce.emplace_back(Force(0, xx));
  // Force( 1/2, xx)
  auto q = xx.linear ? 1.0 : 1.0 + mD_n[0] * mD_n[0];
  double a = mASetForce[0] - 2.0 * xx.xi * mVel_n[0] - mD_n[0] * q;
  mAccel_n.emplace_back(a);
  // the time integrator has to be set up for the periodic loading function.
  double period = 2 * M_PI / xx.circular_frequency;
  double cycles = mTime_integ.TotalDuration() / period;
  mNumCycles = round(cycles);
  assert(std::abs(static_cast<double>(mNumCycles) - cycles) < 1.e-8);
  int n = mTime_integ.NSteps();
  assert(n % mNumCycles == 0);
  return EXIT_SUCCESS;
}

// force for the trapezoid rule,  f_k+1/2
double SpaceTimeTransient::Force(int step, HB x) {
  int s = mTime_integ.StepsTaken();
  assert(s <= step && step <= s + 1);
  int n = static_cast<int>(mASetForce.size());
  if (step == n) {
    double time = mTime_integ.TimeAtStartOfStep();
    time += mTime_integ.mCurrentStepSize;
    auto f = x.amp * sin(x.circular_frequency * time);
    mASetForce.emplace_back(f);
    return f;
  }
  return mASetForce[step];
}

std::vector<double> SpaceTimeTransient::Load(int step, HB x) {
  std::vector<double> rhs = {0.0, Force(step, x)};
  return rhs;
}

// y' = A y,   y' = [u,v]
std::vector<double> SpaceTimeTransient::OdeMatrix(double xi, double u) {
  double q = 1.0 + u * u;
  std::vector<double> A = {0.0, -q, 1.0, -2.0 * xi};
  return A;
}
// rhs = yo + (h/2)(A(yo)yo + f0 + f1)
std::vector<double> SpaceTimeTransient::TrapezoidRhs(HB xx) {
  int step = mTime_integ.StepsTaken();
  std::vector<double> yo = {mD_n[step], mVel_n[step]};  // A_ij y_j
  double u = xx.linear ? 0.0 : mD_n[step];
  std::vector<double> A = OdeMatrix(xx.xi, u);
  std::vector<double> Ay = mat_vec_2x2(yo, A);
  std::vector<double> fo = Load(step, xx);      // mASetForce[0];
  std::vector<double> f1 = Load(step + 1, xx);  // mASetForce[1];
  std::vector<double> rhs = yo;
  double h = mTime_integ.StepSize();
  for (int i = 0; i < 2; i++) {
    rhs[i] += (Ay[i] + fo[i] + f1[i]) * 0.5 * h;
  }
  return rhs;
}

// r(y1) = (I - (h/2) A(y1)) y1 - rhs
std::vector<double> SpaceTimeTransient::NonlinearResidual(HB xx) {
  int step = mTime_integ.StepsTaken();
  int jj = static_cast<int>(mD_n.size());
  assert(step + 1 < jj);
  jj = static_cast<int>(mVel_n.size());
  assert(step + 1 < jj);
  double y1[2] = {mD_n[step + 1], mVel_n[step + 1]};
  double h = mTime_integ.StepSize();
  double u = xx.linear ? 0.0 : mD_n[step + 1];
  std::vector<double> A = OdeMatrix(xx.xi, u);
  for (int i = 0; i < 4; i++) {
    A[i] *= -(0.5 * h);
  }
  A[0] += 1.0;
  A[3] += 1.0;
  std::vector<double> resid = TrapezoidRhs(xx);
  resid[0] = -resid[0];
  resid[1] = -resid[1];
  resid[0] = resid[0] + A[0] * y1[0] + A[2] * y1[1];
  resid[1] = resid[1] + A[1] * y1[0] + A[3] * y1[1];
  return resid;
}

// r(y) = (I - A(y)h/2) y - rhs
//  J = (I- A h/2) - A'(y) y h/2
std::vector<double> SpaceTimeTransient::JacobianMatrix(HB xx) {
  int step = mTime_integ.StepsTaken();
  int jj = static_cast<int>(mD_n.size());
  assert(step + 1 < jj);
  double xi = xx.xi;
  double u = xx.linear ? 0.0 : mD_n[step + 1];
  double r = 1.0 + 3.0 * u * u;
  double h = mTime_integ.StepSize();
  std::vector<double> J = {1.0, r * 0.5 * h, -0.5 * h, 1.0 + xi * h};
  return J;
}

std::vector<double> SpaceTimeTransient::fdJacobianMatrix(HB xx, double tau) {
  assert(!xx.linear);
  std::vector<double> ro = NonlinearResidual(xx);
  int step = mTime_integ.StepsTaken();
  int jj = static_cast<int>(mD_n.size());
  assert(step + 1 < jj);
  double u = mD_n[step + 1];
  jj = static_cast<int>(mVel_n.size());
  assert(step + 1 < jj);
  double v = mVel_n[step + 1];
  mD_n[step + 1] = u + tau;
  std::vector<double> J1 = NonlinearResidual(xx);
  mD_n[step + 1] = u;
  mVel_n[step + 1] = v + tau;
  std::vector<double> J2 = NonlinearResidual(xx);
  mVel_n[step + 1] = v;

  for (int i = 0; i < 2; i++) {
    J1[i] -= ro[i];
    J2[i] -= ro[i];
    J1[i] /= tau;
    J2[i] /= tau;
  }
  std::vector<double> J = {J1[0], J1[1], J2[0], J2[1]};
  return J;
}

int SpaceTimeTransient::TrapezoidPredictor() {
  int step = mTime_integ.StepsTaken();
  double h = mTime_integ.StepSize();
  int sz = static_cast<int>(mD_n.size());
  assert(step + 1 == sz);
  // Temporary hack, use D,V,A, F as space time data.
  mD_n.emplace_back(mD_n[step] + (mVel_n[step] + 0.5 * h * mAccel_n[step]) * h);
  mVel_n.emplace_back(mVel_n[step] + mAccel_n[step] * h);
  return (step + 2 == sz ? EXIT_SUCCESS : EXIT_FAILURE);
}

int SpaceTimeTransient::NonlinearSolveAtEachStep(HB xx) {
  int num_iter = (xx.linear ? 1 : 2);
  for (int iter = 0; iter < num_iter; iter++) {
    std::vector<double> resid = NonlinearResidual(xx);
    std::vector<double> J = JacobianMatrix(xx);
    resid[0] = -resid[0];
    resid[1] = -resid[1];
    std::vector<double> sol = solve_2by2(resid, J);
    std::vector<double> rmJs = residual_2by2(resid, J, sol);
    assert(std::abs(rmJs[0]) + std::abs(rmJs[1]) < 1.e-12);
    int jj = static_cast<int>(mD_n.size());
    int step = mTime_integ.StepsTaken();
    assert(step + 1 < jj);
    mD_n[step + 1] += sol[0];
    mVel_n[step + 1] += sol[1];
  }
  return EXIT_SUCCESS;
}

int SpaceTimeTransient::Solve(HB xx) {
  double h = mTime_integ.StepSize();
  for (int step = 0; step < mTime_integ.NSteps(); step++) {
    // Temporary hack, use D,V,A, F as space time data.
    TrapezoidPredictor();
    NonlinearSolveAtEachStep(xx);

    if (step > 0 && mVel_n[step] >= 0.0 && mVel_n[step + 1] <= 0.0) {
      find_max_disp(step);
    }
    double u = mD_n[step + 1];
    double q = (xx.linear ? 1.0 : 1.0 + u * u);
    double ddot = mASetForce[step + 1] - 2.0 * xx.xi * mVel_n[step + 1] - u * q;
    mAccel_n.emplace_back(ddot);
    int sza = static_cast<int>(mAccel_n.size());
    assert(sza == step + 2);
    mTime_integ.IncrementTimeBeforeCallingTakeNextStep();
    mTime_integ.TakeNextStep();
  }
  return write_trajectory_to_disk(h);
}

//  u'' + xi u' + u = a sin(wt)
int SpaceTimeTransient::Wip(HB xx) {
  double w = xx.circular_frequency;

  // L = w^2 D2 + 2 xi w D + I
  std::vector<double> Lvec = getIdentity(xx.order);
  std::vector<double> D2 = getDsquared(xx.order);
  axpy(D2, w * w, Lvec);
  std::vector<double> D = getD(xx.order);
  axpy(D, 2.0 * xx.xi * w, Lvec);
  auto n = 2 * xx.order + 1;  // mTime_integ.NSteps()
  // print_square_matrix_col(Lvec);
  std::vector<double> rhs(n, 0.0);
  double h = 2.0 * M_PI / static_cast<double>(n);  // omega * dt
  double a = xx.amp;
  for (int i = 0; i < n; i++) {
    rhs[i] = a * sin(static_cast<double>(i) * h);  // mASetForce[i]
  }
  // std::cout << " norm of rhs is " << normInf(rhs) << "\n";
  //  Begin solve  L x = rhs,  where x ~ mD_n
  std::pair<int, int> szA = {n, n};
  std::vector<std::vector<double>> A = matrixOfVector(true, szA, Lvec);
  // std::cout << " norm of L is " << norm1(szA,A) << "\n";
  std::vector<std::vector<double>> Q;  // Q'Q = I
  std::vector<std::vector<double>> R;  // A = QR
  // QRdecomosition is the simple interface to Householder QR.
  // The efficient interface is called HQR2.
  // print_matrix("A ", A);
  QRdecomposition(A, Q, R);
  // print_matrix("orth", Q);
  // print_matrix("triu", R);

  // solve will fail if a diagonal of R is sufficiently small
  auto Qtrhs = gemvt(szA, 1.0, Q, rhs, 0);
  auto sol = apply_inv_triu(szA, R, Qtrhs);
  auto res_in = apply_triu(szA, R, sol);
  axpy(Qtrhs, -1.0, res_in);
  // std::cout << " inner residual norm is " <<
  const double tol = 1.e-14;
  assert(normInf(res_in) < tol);
  auto res = gemv(szA, 1.0, A, sol, 0);
  axpy(rhs, -1.0, res);
  // std::cout << " residual norm is " <<
  assert(normInf(res) < tol);
  // end solve  L x = rhs
  double max = 0.0;
  int argmax = 0;
  for (int i = 0; i < n; i++) {
    if (std::abs(sol[i]) > max) {
      argmax = i;
      max = std::abs(sol[i]);
    }
  }
  // double dt = h/w; // mTime_integ.StepSize();
  std::cout << "diag R >= " << std::abs(min_abs_diag(szA, R)) << "   sol <= sol_" << argmax << " = "
            << max << "\n";
  return EXIT_FAILURE;
}

int SpaceTimeTransient::write_trajectory_to_disk(double step_size) {
  if (!write_to_file) {
    return EXIT_SUCCESS;
  }
  size_t szd = mD_n.size();
  std::vector<double> time(szd, 0.0);
  for (size_t step = 1; step < szd; step++) {
    time[step] = time[step - 1] + step_size;
  }
  Print4Vectors(time, mD_n, mVel_n, mAccel_n);
  return EXIT_SUCCESS;
}

void SpaceTimeTransient::StoreSolution() {
  int restartstep = mTime_integ.StepsTaken() + 1;
  int start = mNeq * restartstep;
  std::copy(mD_n.begin(), mD_n.end(), std::begin(mSpaceTimeDisp) + start);
  std::copy(mVel_n.begin(), mVel_n.end(), std::begin(mSpaceTimeVel) + start);
  std::copy(mAccel_n.begin(), mAccel_n.end(), std::begin(mSpaceTimeAccel) + start);
  std::copy(mASetForce.begin(), mASetForce.end(), std::begin(mSpaceTimeForce) + start);
}

void SpaceTimeTransient::SetInitialCondition(std::vector<double> &deformation,
                                             std::vector<double> &velocity) {
  // std::vector<double> force;
  // std::copy(force.begin(), force.end(), std::begin(mSpaceTimeForce));
  mD_n.emplace_back(deformation[0]);
  mVel_n.emplace_back(velocity[0]);
  std::copy(deformation.begin(), deformation.end(), std::begin(mSpaceTimeDisp));
  std::copy(velocity.begin(), velocity.end(), std::begin(mSpaceTimeVel));
}

void SpaceTimeTransient::PrintMe() { mTime_integ.PrintMe(); }
void SpaceTimeTransient::OutputTrajectory(bool default_is_false) {
  write_to_file = default_is_false;
}

void SpaceTimeTransient::find_max_disp(int step) {
  double um = mD_n[step];
  double up = mD_n[step + 1];
  double vm = mVel_n[step];
  double vp = mVel_n[step + 1];
  double dvm = mAccel_n[step];
  double dvp = mAccel_n[step + 1];
  double h = mTime_integ.StepSize();
  auto tu = init_disp(um, up, vm, vp, dvm, dvp, h);
  int steps_per_cycle = mTime_integ.NSteps() / mNumCycles;
  int n = mTime_integ.StepsTaken() % steps_per_cycle;
  double time = (static_cast<double>(n) + 0.5) * h + tu.first;
  std::cout << "time " << time << " displacement " << tu.second << "\n";
}
void SpaceTimeTransient::Print4Vectors(const std::vector<double> &t,
                                       const std::vector<double> &d,
                                       const std::vector<double> &v,
                                       const std::vector<double> &a) {
  std::ofstream fout;
  fout.open("dva.dat");
  size_t n = d.size();
  for (size_t i = 0; i < n; i++) {
    // fout << std::setw(23) << std::setprecision(15);
    fout << t[i] << "     " << d[i] << "  " << v[i] << "  " << a[i] << "\n";
  }
  fout.close();
}
