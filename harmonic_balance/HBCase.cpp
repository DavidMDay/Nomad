#include <array>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>  // inner product ?
#include <numeric>
#include "FourierCoefficients.hpp"
#include "TemporalParameters.hpp"
#include "FixedStepIntegrator.hpp"
#include "Tests.hpp"
#include "SpaceTimeTransient.hpp"

TemporalParameters get_my_tp(){
  TemporalParameters t;
  t.mTimeStep.emplace_back(0.1);
  t.mNumSkip.emplace_back(1);
  t.mNumStep.emplace_back(2);
  return t;;
}


TemporalParameters get_a_different_tp(){  // unsued
  TemporalParameters t = get_my_tp();
  t.mTimeStep[0]=0.025;
  t.mNumStep[0]=10;
  return t;
}


material get_my_material()  // unused
{
  material x;
  x.E = 2.7e7;
  x.density = 7.5e-4;
  x.poissons_ratio = 0.3;
  return x;
}

HB get_another_hb() { // unused
  HB xx;
  xx.circular_frequency = 25.132741228718346;
  return xx;
}

int TestBoundarycondition() {
   // const salinas_sol_case* IVP = SalinasGlobals::singleton()->Multicase()->ThisSolution();
   // const harmonicbalance_sol* ode = dynamic_cast<const harmonicbalance_sol*>(IVP);
    double dt = ode->flavor.time_step;
    int nstep = ode->flavor.nsteps;
    double period = 2.0*M_PI/ode->flavor.circular_frequency;
    double duration = dt * static_cast<double>(nstep);
    EXPECT_NEAR(period, duration, 1.0e-9);

    std::valarray<double> force(nstep);
    double time = 0.0;
    for (int step = 0; step < nstep; step++) {
      double appliedForce = 1.0e2 * sin(8.0 * M_PI * time);
      force[step] = source->Value(time);
      EXPECT_NEAR(appliedForce, force[step], 1.0e-9);
      time += dt;
    }
    LOADREGION region = load->RegionType();
    assert(region == LOADREGION::NODE_SET);
    int four = 4; // a quadrilateral
    assert(4 == four);
    int nodeoffset = set->NodeId(onlyOne);
    assert(24 == nodeoffset);  // 6 elements

    std::valarray<std::complex<double>> phat(0.0,nstep);
    FourierCoefficients::Dft(force, phat);

    double au = 50.0 * static_cast<double>(nstep);
    EXPECT_NEAR(-au, imag(phat[1]), 1.0e-9);
    EXPECT_NEAR( au, imag(phat[nstep-1]), 1.0e-9);
}

int TestSolve() {
    IVP->Initialize();
    harmonicbalance_sol* unit_testable = dynamic_cast<harmonicbalance_sol*>(IVP);

    std::pair<mtk::Vec3<double>, mtk::Vec3<double>> box =  unit_testable->BoundingBox();
    std::cout << " bounding box"<< box.second << "\n";
    double origin_size = mtk::infinity_norm( box.first );
    EXPECT_NEAR(0.0, origin_size, 1.0e-8);
    mtk::Vec3<double> point = {1.0,1.0,6.0};
    double error = mtk::infinity_norm( point - box.second );
    EXPECT_NEAR(0.0, error, 1.0e-8);

    unit_testable->Solve(sal_globals->Multicase());
    std::stringstream hay;
    hay << "solve once and print 1-norm of Fourier coefficient at each dof\n";
    const FourierCoefficients& forceHat = unit_testable->GetLoadHat();
    const FourierCoefficients& dispHat = unit_testable->GetSolHat();
    FourierCoefficients forceHatNorm = forceHat.norm1();
    FourierCoefficients dispHatNorm = dispHat.norm1();
    forceHatNorm.Print(0,hay);
   // dispHatNorm.Print(0,hay);
   std::cout << hay.str();
   
    std::vector<double> f = forceHat.SineField(1);
    std::vector<double> u = dispHat.SineField(1);
    std::vector<double> Ku = unit_testable->Apply( Kaa, u.data());
    std::vector<double> Mu = unit_testable->Apply( Maa, u.data());

    std::vector<double> r(f);
    int neq = dbm->GetValue<int>("NEQ");
    for( int i=0; i<neq; i++ ) {
        r[i] -= Ku[i];
    }
    double f2 = std::inner_product(f.begin(), f.end(), f.begin(), 0.0); // 2e4
    double r2 = std::inner_product(r.begin(), r.end(), r.begin(), 0.0); // 1e-11
    assert(r2 < f2);
    double m2 = std::inner_product(Mu.begin(), Mu.end(), Mu.begin(), 0.0);// 2.5 e-18
    double rm = std::inner_product(r.begin(), r.end(), Mu.begin(), 0.0); // -5 e-15
    // r = f - (K - sM)u =  f - Ku + theta (Mu)
    // The Fourier series is at 0, theta, 2 theta, ....  .   
    double theta = -rm/m2;   //  5.0e-15 /2.5 e-18  = 2e3

    std::cout << f2 << "  " << r2 << "  " << m2 << "  " << rm << "   " << theta << "\n";

    IVP->Finalize();
}

int TestReanimated) {
  SalinasGlobals* sal_globals = SalinasGlobals::singleton();
  salinas_sol_case* IVP = SalinasGlobals::singleton()->Multicase()->ThisNonConstSolution();
  IVP->Initialize(sal_globals->Multicase());
  harmonicbalance_sol* unit_testable = dynamic_cast<harmonicbalance_sol*>(IVP);
  unit_testable->Solve(sal_globals->Multicase());
  const SpaceTimeTransient* stt = unit_testable->GetTransDriver();
  salinas_dbmgr* dbm = SalinasGlobals::singleton()->Multicase()->Dbm();
  int neq = dbm->GetValue<int>("NEQ");
  int numLoadedNodes = 4;
  double time{0.0};
  const double dt{0.025};
  const int nstep{10};
  for (int step = 0; step < nstep; step++) {
      double loadSum = 0.0;
      for (int idof = step * neq; idof < (step + 1) * neq; ++idof) {
        loadSum += stt->mSpaceTimeForce[idof];
      }
      double appliedForce = numLoadedNodes * 1.0e2 * sin(8.0 * M_PI * time);
      EXPECT_NEAR(appliedForce, loadSum, 1.0e-9);
      time += dt;
  }

  IVP->Finalize();
}
