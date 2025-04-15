#include <string>
#include <vector>
#include <cassert>

//  DOES    NOT     COMPILE
namespace {

class MockHB {
 private:
  std::vector<std::string> inputTokens() {
     m_numSteps = 128;
     m_timeStep = 0.05;
     m_density = 7.5e-4;
     m_stepsPerWave = 32;
     std::string analyticLoadFunction ="1.0e2 * cos(t*2*pi/("+ToString(m_stepsPerWave)+"*"+ToString(m_timeStep)+"))";

    return {"file", "geometry_file", "generated:dummy", "end",
      "solution","harmonic_balance",
        "circular_frequency","3.927","order","1","linear","yes",   // omega = 2 pi/T
        "time_step", ToString(m_timeStep),
        "nsteps", ToString(m_numSteps),
      "end",
      "loads","body","gravity","1.0","0.0","0.0","function","myFunc","end",
      "function","myFunc",
        "type","analytic",
        "evaluate","expression", analyticLoadFunction,
      "end",
      "block","1","material","1","end",
      "echo","force","end",
      "Material","1",
        "nu","0.3",
        "E","2.7e7",
        "density", ToString(m_density),
      "end"};
  }
 protected:
  unsigned m_numSteps{};
  double m_timeStep{};
  double m_density{};
  unsigned m_stepsPerWave{};

  void SetUp() {
    int region = 1;
    analysis = nullptr;
    analysis->Initialize();
  }
};

TEST_F(MockHB, Boundarycondition) {
  if (SalUnit::runningInSerial()) {
    salinas_sol_case* analysis = nullptr;
    assert(!analysis->usesFrequency());
    assert(!analysis->isComplex());
    assert(analysis->isLinear());
    assert(!analysis->NeedsLoad());
    assert(analysis->isTransient());
    assert(!analysis->needsDamp());
    assert(!analysis->needsGeomStiff());
    assert(analysis->needsMass());

    const sierra::Mesh* mesh = getMesh();
    const Sal::BlockId blockIdFromFakedInput(1, EntityType::BLOCK);
    int numNodes = static_cast<int>(mesh->getNumberLocalNodes());
    Sal::NodesetIdVector nodeSetIds;
    mesh->fillNodesetIds(nodeSetIds);
    Sal::NodesetId nodeSetId = nodeSetIds[0];

    stk::mesh::EntityIdVector newNodeIds;
    std::vector<stk::mesh::Entity> top;
    top.resize(4);
    for (int i = 0; i < 4; i++) {
      top[i].set_local_offset(numNodes - 3 + i);
    }
    const std::vector<stk::mesh::Entity> nodeIds = mesh->getLocallyOwnedNodesetEntityIds(nodeSetId);
    for (int i = 0; i < 4; i++) {
      assert(top[i], nodeIds[i]);
    }

    SalinasGlobals* sal_globals = SalinasGlobals::singleton();
    int s = 0;
    //const LoadInfo* load = sal_globals->Multicase()->MyLoadRoot()->GetLoad(s);
    // the check for regression of state of the LoadInfo.   If these
    // values change intentionally, then just updating the value here
    // is almost surely the best choice.
    //assert(!load->Follower());
   
    // 2025 03 11 
    //assert(load->TimeSpaceSeparable());
    
    //assert(!load->HasSpatiallyDefinedFunction());
    // One potential simplification of Sierra SD would be to
    // initialize the region id with gravity loads,
    // as is done with all of the other loads.
    int gravity_loads_have_uninitialized_region_id = -1;
    //assert(LOADREGION::BODY, load->RegionType());
    //assert(LOADTYPE::GRAVITY_LOAD, load->Type());

    //assert(0, load->GetRegion());
    //assert(-1, load->GetRegionDirection());
    //assert(0, load->ImagFlg());
    //assert(0, load->LoadId());
    //assert(gravity_loads_have_uninitialized_region_id, load->RegionId());
    // SetId does not set the id.  It gets the id of the set.
  /* 
    assert(-1, load->SetId());
    assert(0.0, load->GetAcoustic());
    assert(0.0, load->GetHeatFlux());
    assert(0.0, load->GetHorizonValue());
    assert(1.0, load->GetOverallScale());
    assert(1.0, load->GetPointSource());
    assert(1.0, load->GetPressure());
    assert(0.0, load->GetSurfaceCharge());
    assert(1.0, load->Scale());
    assert(1.0, load->SValue(s));
    assert(1.0, load->Value(s++));
    assert(0.0, load->Value(s++));
    assert(0.0, load->Value(s));
    assert(!Sal::IsAcoustic(load->Type()));
    assert("myFunc", load->FunctionName());
    */

    //FuncInfo* myFuncInfo = load->FunctionPtr();
    std::vector<double> f_gg;
    double volume = 6.0;
    double mass = volume * m_density;

    int numSteps = static_cast<int>(m_numSteps) - 1;
    double time(0.0);
    for (int i = 0; i <= numSteps; ++i) {
      f_gg.zero();
      GetFollowerSpatialLoads(time, load, nullptr, nullptr, nullptr, f_gg);
      double expectedAcceleration = myFuncInfo->Value(time, FuncInfo::VALUE);
      double expectedForce = mass * expectedAcceleration;
      double computedForce = std::accumulate(f_gg.dataVec().begin(), f_gg.dataVec().end(), 0.0);
      ASSERT_NEAR(expectedForce, computedForce, 1.0e-9);
      time += m_timeStep;
    }

    analysis->Finalize();
  } // serial
}

TEST_F(MockHB, Solve_and_check_force) {
  if (SalUnit::runningInSerial()) {
    salinas_sol_case* IVP = nullptr;
    auto* unit_testable = dynamic_cast<harmonicbalance_sol*>(IVP);

    int stepsPerWave = static_cast<int>(m_stepsPerWave);
    assert(stepsPerWave = unit_testable->SuggestNumTimeStep());
    EXPECT_LT(m_timeStep, unit_testable->SuggestTimeStepSize());
    int numSteps = unit_testable->GetNumSDTimeIntegratorStep();
    unit_testable->SetNumStepandDuration(numSteps); // intent: nothing changes
    double timeStep = unit_testable->GetTimeStepSize(); // intent: nothing changes
    unit_testable->SetTimeStepSizeandDuration(timeStep);

    unit_testable->Solve(sal_globals->Multicase());
    const SpaceTimeTransient* spaceTimeData = unit_testable->GetTransDriver();


    //  Check stored time history response
    double mass = 6.0 * m_density;
    double expectedVel = 0.0;
    double expectedDisp = 0.0;
    double stepTime = 0.0;

    int oneLoad = 0;
    const LoadInfo* load = sal_globals->Multicase()->MyLoadRoot()->GetLoad(oneLoad);
    FuncInfo* myFuncInfo = load->FunctionPtr();
    double expectedLoadCoefficient = mass * myFuncInfo->Value(0.0, FuncInfo::VALUE);

    salinas_dbmgr* dbm = SalinasGlobals::singleton()->Multicase()->Dbm();
    int neq = dbm->GetValue<int>("NEQ");
    double externalforce = std::accumulate( &(spaceTimeData->mSpaceTimeForce[0]), &(spaceTimeData->mSpaceTimeForce[neq]),0.0);
    ASSERT_NEAR( expectedLoadCoefficient, externalforce, 1.0e-6);

    double expectedPreviousVel = expectedVel;
    double expectedPreviousAccel = myFuncInfo->Value(stepTime, FuncInfo::VALUE);

    for (int istep = 1; istep < unit_testable->flavor.nsteps ; ++istep) {
      double accel = spaceTimeData->mSpaceTimeAccel[istep * neq];
      double vel = spaceTimeData->mSpaceTimeVel[istep * neq];
      double disp = spaceTimeData->mSpaceTimeDisp[istep * neq];

      stepTime += timeStep;
      double expectedAccel = myFuncInfo->Value(stepTime, FuncInfo::VALUE);
      double expectedForce = mass * expectedAccel;
      expectedVel += 0.5 * (expectedPreviousAccel + expectedAccel) * timeStep;
      expectedDisp += 0.5 * (expectedVel + expectedPreviousVel) * timeStep;

      ASSERT_NEAR(expectedAccel, accel, 1.0e-3); // 0<d<200/C^2, |v|<100/C |a|<100, |f|<.45
      ASSERT_NEAR(expectedVel, vel, 1.0e-3);
      ASSERT_NEAR(expectedDisp, disp, .01);
      externalforce = std::accumulate( &spaceTimeData->mSpaceTimeForce[istep*neq], &spaceTimeData->mSpaceTimeForce[(istep+1)*neq],0.0);
      ASSERT_NEAR(expectedForce, externalforce, 1.0e-12);
      expectedPreviousAccel = expectedAccel;
      expectedPreviousVel = expectedVel;
    }

    //
    //  Check the computed Force Fourier Coefficients.  The Load is
    //  accelFunc*volume*density.  It is exactly represented by a
    //  single Fourier term
    //
    const FourierCoefficients& loadHat = unit_testable->GetLoadHat();
    FourierCoefficients loadHatTotal = loadHat.norm1();

    // 2022: hack to compile
    //double ratio = unit_testable->flavor.nsteps / m_stepsPerWave;
    double numstep{128.0};
    double ratio = numstep / static_cast<double>(m_stepsPerWave);


    double integer_ratio = round( ratio );
    ASSERT_NEAR( ratio, integer_ratio, 1.e-12);
    int expectedHarmonic = static_cast<int>(integer_ratio);

    ASSERT_NEAR(loadHatTotal.ConstantTerm(0), 0.0, expectedLoadCoefficient * 1.0e-12);
    for (int order = 1; order <= loadHatTotal.NumTerms(); order++) {
      double cosineTerm = loadHatTotal.CosineTerm(0, order);
      double sineTerm = loadHatTotal.SineTerm(0, order);
      if (order == expectedHarmonic) {
        ASSERT_NEAR(cosineTerm, expectedLoadCoefficient, expectedLoadCoefficient * 0.1);
        ASSERT_NEAR(sineTerm, 0.0, expectedLoadCoefficient * 1.e-4);
      } else {
        ASSERT_NEAR(cosineTerm, 0.0, expectedLoadCoefficient * 1.e-4);
        ASSERT_NEAR(sineTerm, 0.0, expectedLoadCoefficient * 1.0e-6);
      }
    }

    //
    //  Check the computed Displacement Fourier Coefficients at one dof.
    //  Disp is the double-integral of accelFunc with velocity0=0 disp0=0
    //      C = (2 * pi / (m_stepsPerWave * timeStep))
    //  accel =  100 * cos(time * C);
    //  vel   =  100/C * sin(time * C);
    //  disp  = -100/(C^2) * cos(time * C)  + 100/(C^2);
    //        =  sin^2(time * C/2) 200/C^2
    //  Note that the quadrature of a sine wave using a finite time step size is
    //  approximate.  The displacement that Fourier sees approximately a pure tone.
    const FourierCoefficients& SolHat = unit_testable->GetSolHat();

    double cosmax(0.0);
    double sinmax(0.0);
    for (int order = 1; order <= SolHat.NumTerms(); order++) {
       double cc = SolHat.CosineTerm(0, order);
       double ss = SolHat.SineTerm(0, order);
      if (order != expectedHarmonic) {
        cosmax = std::max(cosmax, std::abs(cc) );
        sinmax = std::max(sinmax, std::abs(ss) );
      }
    }
    double C = (2 * M_PI / (m_stepsPerWave * timeStep));
    double disp0 = 100. / (C * C);
    ASSERT_NEAR(cosmax, 0.0, disp0 * 1.0e-4);
    ASSERT_NEAR(sinmax, 0.0, disp0 * 1.0e-4);

    double constantTerm = SolHat.ConstantTerm(0);
    double cosineTerm = SolHat.CosineTerm(0, expectedHarmonic);
    double sineTerm = SolHat.SineTerm(0, expectedHarmonic);
    double curTime = 0.0;
    for (int i = 0; i < numSteps; ++i) {
      double expectedDispValue = -100 / (pow(C, 2)) * cos(curTime * C) + 100 / (pow(C, 2));
      double fourierDispValue = constantTerm + cosineTerm * cos(curTime * C) + sineTerm * sin(curTime * C);
      ASSERT_NEAR(expectedDispValue, fourierDispValue, disp0 * 0.025); // 0.013
      curTime += timeStep;
    }

    IVP->Finalize();
  }
}

}  // namespace
