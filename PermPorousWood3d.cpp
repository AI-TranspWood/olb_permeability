/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library based on the resolved Rock example:
 *  Copyright (C) 2024 Jan E. Marquardt, Mathias J. Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 *
 *************************************************************************************
 * VTT Technical Research Centre of Finland, 2025
 * Timo Korhonen, olb-1.7
 * Antti Puisto: Port to olb-1.8.1, generalizations of the code, boundary implementation
 * Project: EU-AI-Transpwood
 * 
 * Use the AITW python generated birch microstructure to construct the geometry
 * Use the AITW filter script to convert the image stack to vti.
 *
 * Permeability calculation (similar to porous Rock case example, exposed viscosity to allow
 * simulations of effective permeability)
 *
 * This folder/file: Runs OpenLB simulation to calculate the permeability in the selected direction
 *
 *
 *
 *************************************************************************************
 *
 *
 * This uses GuoZhaoBGKdynamics for whole domain exept the BoundingBox boundaries.
 * The wood structure is read from the VTI file and its porosity is set to <dSolid>.
 * The permeability is calculated using two methods, one using flow measurement,
 * the other is using overall velocity average. These are quite close to each other,
 * the flow method should be more correct as it takes the flow at a slice going through
 * the wood sample only.
 *
 * Model: Treat the wood cell structure as OpenLB porous model, so
 *        that monomer can avoind "dead ends" of the cells and also can
 *        move other directions than the available "pipes".   *
*/

#include <olb.h>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T          = FLOATING_POINT_TYPE;   // Provided by OpenLB (float/double)
// PorousBGKdynamics/POROSITY dropped in favor of GuoZhaoBGKdynamics, whose Darcy
// permeability (K) is set directly (Kin/h^2) rather than via a (1-d) transform
// with an h^2*nu*tau floor. 
using DESCRIPTOR = D3Q27<descriptors::FORCE, descriptors::EPSILON, descriptors::K,
                          descriptors::NU, descriptors::BODY_FORCE>;

// -----------------------------------------------------------------------------
// Data containers
// -----------------------------------------------------------------------------

// Configuration and inputs parsed from CLI
struct SimulationConfig {
  // I/O
  std::string vtiFile;
  std::string arrayName;

  // Physical & numerical controls
  int   resolution           = 0;       // lattice resolution N
  T     tau                  = T{0};    // relaxation time
  T     InletPressure         = T{0};    // Pa per meter (dp / L)
  T     scalingFactor        = T{1};    // VTI geometry scale
  T     wallPermeability     = T{1e-16}; // target physical Darcy permeability of the
                                          // wood material (material 5), in m^2. Set
                                          // directly 
  T     wallEpsilon          = T{1.0};   // void fraction used in the GuoZhao force term
                                          // for material 5. 1.0 = pure drag term, no
                                          // volume-averaging correction (matches the
                                          // Spaid-Phelan-equivalent configuration).
  T     kinematicViscosity   = T{1e-6};
  T     fluidDensity         = T{1000};
  T     uPhys                = T{0};    // physical velocity magnitude at outlet (m/s)
  T     criterium             = T{1e-5};
  int   flowDirection        = 0;       // 0=x, 1=y, 2=z
  int   wallType              = 0;       // 0: periodic; 1: free-slip; 2: bounce-back
  bool  uniformGuoZhao        = false;   // DIAGNOSTIC: if true, use GuoZhaoBGKdynamics

  // Helper velocity vector (aligned with flowDirection)
  Vector<T,3> uVector{T{0}, T{0}, T{0}};

// Derived convenience
  T dynamicViscosity() const { return kinematicViscosity * fluidDensity; }
};

// Geometry derived from VTI (and padded in flow direction for BCs)
struct GeometryData {
  Vector<T,3> extent{};     // bbox with padding (main domain)
  Vector<T,3> origin{};     // (kept for completeness)
  Vector<T,3> extentVTI{};  // scaled VTI extent (with padding)
  Vector<T,3> originVTI{};
};

// Values produced during the simulation (for logging/consumers)
struct SimulationResults {
  Vector<T,3> averageVelocity{T{0},T{0},T{0}};
  T           permeability = T{0};
};

// Long-lived objects that must persist
struct RuntimeContext {
  std::shared_ptr<BlockVTIreader3D<T,T>> vtiReader;
};

// -----------------------------------------------------------------------------
// Functions
// -----------------------------------------------------------------------------

std::shared_ptr<IndicatorBlockData3D<T>>
generateIndicatorFromVTI(const SimulationConfig& cfg,
                         int flowDirection,
                         GeometryData& geom,
                         RuntimeContext& ctx)
{
  const T sourceScale = cfg.scalingFactor;

  ctx.vtiReader = std::make_shared<BlockVTIreader3D<T,T>>(cfg.vtiFile, cfg.arrayName);

  auto           cuboidSample     = ctx.vtiReader->getCuboid();
  T              deltaRsample     = cuboidSample.getDeltaR() * sourceScale;
  Vector<int,3>  extentSample     = cuboidSample.getExtent();
  Vector<T,3>    originSamplePhys = cuboidSample.getOrigin() * sourceScale;
  Vector<T,3>    extentSamplePhys = { deltaRsample * T(extentSample[0]),
                                      deltaRsample * T(extentSample[1]),
                                      deltaRsample * T(extentSample[2]) };

  // Pad in the flow direction to provide buffer for inlet/outlet
  extentSamplePhys[flowDirection] += 10 * deltaRsample;
  originSamplePhys[flowDirection] -= 5 * deltaRsample;

  // Publish geometry
  geom.extentVTI = extentSamplePhys;
  geom.originVTI = originSamplePhys;
  geom.extent    = extentSamplePhys;
  geom.origin    = originSamplePhys;

  return std::make_shared<IndicatorBlockData3D<T>>(
      ctx.vtiReader->getBlockData(), extentSamplePhys, originSamplePhys, deltaRsample, true);
}

void prepareGeometry(const SimulationConfig& cfg,
                     const UnitConverter<T,DESCRIPTOR>& converter,
                     IndicatorF3D<T>& indicator,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // 1) Start as solid
  superGeometry.rename(0, 5);              // ALL interior -> Material 5
  // 2) Where indicator TRUE -> fluid
  superGeometry.rename(5, 1, indicator);   // Material 5 -> 1 where fluid

  // 3) Inlet/outlet from fluid bounds
  Vector<T,3> minPos = superGeometry.getStatistics().getMinPhysR(1);
  Vector<T,3> maxPos = superGeometry.getStatistics().getMaxPhysR(1);

  // Inflow (Material 3)
  Vector<T,3> inflowOrigin = minPos;
  Vector<T,3> inflowExtent = maxPos - minPos;
  inflowOrigin[cfg.flowDirection] = minPos[cfg.flowDirection] - T{1.5} * converter.getPhysDeltaX();
  inflowExtent[cfg.flowDirection] = T{1.0} * converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(inflowExtent, inflowOrigin);
  superGeometry.rename(1, 3, inflow);
  superGeometry.rename(5, 3, inflow);

  // Outflow (Material 4)
  Vector<T,3> outflowOrigin = minPos;
  Vector<T,3> outflowExtent = maxPos - minPos;
  outflowOrigin[cfg.flowDirection] = maxPos[cfg.flowDirection] - T{0.5} * converter.getPhysDeltaX();
  outflowExtent[cfg.flowDirection] = T{0.75} * converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(outflowExtent, outflowOrigin);
  superGeometry.rename(1, 4, outflow);
  superGeometry.rename(5, 4, outflow);

  // Clean inner materials
  superGeometry.innerClean();

  // Stats
  clout << std::endl;
  clout << "========== Final Material Distribution ==========" << std::endl;
  superGeometry.print();
  clout << "Material 0 (exterior):   " << superGeometry.getStatistics().getNvoxel(0) << std::endl;
  clout << "Material 1 (fluid):      " << superGeometry.getStatistics().getNvoxel(1) << std::endl;
  clout << "Material 2 (walls):      " << superGeometry.getStatistics().getNvoxel(2) << std::endl;
  clout << "Material 3 (inflow):     " << superGeometry.getStatistics().getNvoxel(3) << std::endl;
  clout << "Material 4 (outflow):    " << superGeometry.getStatistics().getNvoxel(4) << std::endl;
  clout << "Material 5 (porous):     " << superGeometry.getStatistics().getNvoxel(5) << std::endl;
  clout << "=================================================" << std::endl;

  const char* wallTypeStr[] = {"PERIODIC (no walls)", "FREE-SLIP walls", "BOUNCE-BACK walls"};
  clout << "Wall type: " << wallTypeStr[cfg.wallType] << std::endl;
  clout << "Flow direction: "
        << (cfg.flowDirection == 0 ? "X" : cfg.flowDirection == 1 ? "Y" : "Z") << std::endl;

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(const SimulationConfig& cfg,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    IndicatorF3D<T>& indicator,
                    SuperLattice<T,DESCRIPTOR>& sLattice,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  const T h_SP  = converter.getPhysDeltaX();
  const T nu_SP = (cfg.tau - T{0.5}) / T{3};   // lattice kinematic viscosity

  if (cfg.uniformGuoZhao) {
    clout << "*** DIAGNOSTIC MODE: uniform GuoZhaoBGKdynamics (fluid + wood) ***" << std::endl;

    auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
    sLattice.defineDynamics<GuoZhaoBGKdynamics<T,DESCRIPTOR>>(bulkIndicator);
    clout << "GuoZhaoBGKdynamics applied to Materials 1, 3, 4, 5" << std::endl;

    // Walls if enabled
    if (cfg.wallType == 1) {
      clout << "Setting FREE-SLIP boundaries on Material 2..." << std::endl;
      boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
    } else if (cfg.wallType == 2) {
      clout << "Setting BOUNCE-BACK boundaries on Material 2..." << std::endl;
      boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
    } else {
      clout << "No wall boundaries (periodic BC)" << std::endl;
    }

    const T fluidK = T{1.0};   // m^2, deliberately huge -> ~negligible drag (nu/K ~ 0)

    clout << "\n========== Setting GuoZhao Fields (all materials) ==========\n";
    clout << "  Materials 1,3,4 (fluid): epsilon=1.0, K=" << fluidK << " m^2 (negligible drag)\n";
    clout << "  Material 5 (wood):       epsilon=" << cfg.wallEpsilon
          << ", K=" << cfg.wallPermeability << " m^2" << std::endl;

    AnalyticalConst3D<T,T> fluidEpsilon(T{1.0});
    AnalyticalConst3D<T,T> fluidNu(nu_SP);
    AnalyticalConst3D<T,T> fluidKfield(fluidK / (h_SP * h_SP));
    AnalyticalConst3D<T,T> fluidBodyForce(Vector<T,3>{T{0},T{0},T{0}});

    for (int mat : {1,3,4}) {
      sLattice.defineField<descriptors::EPSILON>(superGeometry, mat, fluidEpsilon);
      sLattice.defineField<descriptors::NU>(superGeometry, mat, fluidNu);
      sLattice.defineField<descriptors::K>(superGeometry, mat, fluidKfield);
      sLattice.defineField<descriptors::BODY_FORCE>(superGeometry, mat, fluidBodyForce);
    }

    AnalyticalConst3D<T,T> woodEpsilon(cfg.wallEpsilon);
    AnalyticalConst3D<T,T> woodNu(nu_SP);
    AnalyticalConst3D<T,T> woodK(cfg.wallPermeability / (h_SP * h_SP));
    AnalyticalConst3D<T,T> woodBodyForce(Vector<T,3>{T{0},T{0},T{0}});

    sLattice.defineField<descriptors::EPSILON>(superGeometry, 5, woodEpsilon);
    sLattice.defineField<descriptors::NU>(superGeometry, 5, woodNu);
    sLattice.defineField<descriptors::K>(superGeometry, 5, woodK);
    sLattice.defineField<descriptors::BODY_FORCE>(superGeometry, 5, woodBodyForce);

    clout << "==============================================================\n";

  } else {
    // DEFAULT (production) MODE: split dynamics -- plain BGKdynamics for fluid,
    // GuoZhaoBGKdynamics only for the wood material.
    auto fluidIndicator = superGeometry.getMaterialIndicator({1,3,4});
    sLattice.defineDynamics<BGKdynamics<T,DESCRIPTOR>>(fluidIndicator);
    clout << "BGKdynamics applied to Materials 1, 3, 4 (fluid)" << std::endl;

    auto woodIndicator = superGeometry.getMaterialIndicator({5});
    sLattice.defineDynamics<GuoZhaoBGKdynamics<T,DESCRIPTOR>>(woodIndicator);
    clout << "GuoZhaoBGKdynamics applied to Material 5 (wood)" << std::endl;

    // Walls if enabled
    if (cfg.wallType == 1) {
      clout << "Setting FREE-SLIP boundaries on Material 2..." << std::endl;
      // Fallback to BounceBack if dedicated free-slip is not available
      boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
    } else if (cfg.wallType == 2) {
      clout << "Setting BOUNCE-BACK boundaries on Material 2..." << std::endl;
      boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
    } else {
      clout << "No wall boundaries (periodic BC)" << std::endl;
    }

    // GuoZhao permeability fields -- material 5 (wood) only. Fluid materials
    // (1,3,4) never use these, since BGKdynamics doesn't read them.
    clout << "\n========== Setting Wood Permeability Fields (Material 5) ==========\n";
    clout << "  wallPermeability (Kin): " << cfg.wallPermeability << " m^2\n";
    clout << "  wallEpsilon:            " << cfg.wallEpsilon << std::endl;

    AnalyticalConst3D<T,T> woodEpsilon(cfg.wallEpsilon);
    AnalyticalConst3D<T,T> woodNu(nu_SP);
    AnalyticalConst3D<T,T> woodK(cfg.wallPermeability / (h_SP * h_SP));  // K in lattice units
    AnalyticalConst3D<T,T> woodBodyForce(Vector<T,3>{T{0},T{0},T{0}});

    sLattice.defineField<descriptors::EPSILON>(superGeometry, 5, woodEpsilon);
    sLattice.defineField<descriptors::NU>(superGeometry, 5, woodNu);
    sLattice.defineField<descriptors::K>(superGeometry, 5, woodK);
    sLattice.defineField<descriptors::BODY_FORCE>(superGeometry, 5, woodBodyForce);
    // Note: descriptors::FORCE is recomputed from (u, epsilon, k, nu, bodyF) every
    // collision step by GuoZhaoForcing, so it needs no initial value here.

    clout << "  Wood nodes (material 5): " << superGeometry.getStatistics().getNvoxel(5) << std::endl;
    clout << "====================================================================\n";
  }

  // Inflow: LocalPressure; Outflow: LocalVelocity. For the most stable boundary conditions
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 4);
  clout << "BCs: Material 3 = LocalPressure (inflow), Material 4 = LocalVelocity (outflow)\n";

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF(T{1});
  AnalyticalConst3D<T,T> uF(Vector<T,3>{T{0},T{0},T{0}});

  if (cfg.wallType == 0) {
    sLattice.defineRhoU(superGeometry.getMaterialIndicator({1,3,4,5}), rhoF, uF);
    for (int i: {1,3,4,5}) {
      sLattice.iniEquilibrium(superGeometry, i, rhoF, uF);
    }
  } else {
    sLattice.defineRhoU(superGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
    for (int i: {1,2,3,4,5}) {
      sLattice.iniEquilibrium(superGeometry, i, rhoF, uF);
    }
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Smoothly ramp inlet pressure; set outlet velocity to user-specified u
void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& sLattice,
                       const UnitConverter<T,DESCRIPTOR>& converter,
                       const SimulationConfig& cfg,
                       const GeometryData& geom,
                       int iT, int iTmaxStart,
                       SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  // --- Outlet velocity (constant) in lattice units ---
  Vector<T,3> uPhysVec{T{0},T{0},T{0}};
  uPhysVec[cfg.flowDirection] = cfg.uPhys;

  // Convert physical velocity magnitude to lattice units
  const T uLatMag = converter.getLatticeVelocity(cfg.uPhys); 
  Vector<T,3> uLatVec{T{0},T{0},T{0}};
  uLatVec[cfg.flowDirection] = uLatMag;

  AnalyticalConst3D<T,T> uOut(uLatVec);
  sLattice.defineU(superGeometry, 4, uOut); // Material 4 = LocalVelocity

  // --- Inlet pressure ramp (to prevent start-up shock) ---
  if (iT <= iTmaxStart) {
    const int iTupdate = std::max(iTmaxStart / 200, 1);
    if (iT % iTupdate == 0) {
      PolynomialStartScale<T,int> StartScale(iTmaxStart, T(1));
      int iTvec[1] = {iT}; T frac[1] {};
      StartScale(frac, iTvec);

      const T currentPressure = frac[0] * cfg.InletPressure * geom.extentVTI[cfg.flowDirection];

      AnalyticalConst3D<T,T> rhoIn(converter.getLatticeDensityFromPhysPressure(currentPressure));
      sLattice.defineRho(superGeometry, 3, rhoIn); // Inlet only

      clout << "step= " << iT << "; ramped P_inlet= " << currentPressure << " Pa" << std::endl;

      sLattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(ProcessingContext::Simulation);
    }
  }
}

// Read the outputs
void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                const UnitConverter<T,DESCRIPTOR>& converter,
                const SimulationConfig& cfg,
                const GeometryData& geom,
                SimulationResults& results,
                int iT, T maxPhysT,
                util::Timer<T>& timer,
                bool converged,
                int iTmaxStart,
                SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "getResults");

  // Writers & functors
  SuperVTMwriter3D<T>                       vtmWriter("woodstructure");
  SuperLatticePhysVelocity3D<T,DESCRIPTOR>  velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T,DESCRIPTOR>  pressure(sLattice, converter);

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  const int statIter = converter.getLatticeTime(maxPhysT * T{0.01});
  const int vtkIter  = converter.getLatticeTime(maxPhysT * T{0.05});

  if (iT == 0) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR>   rank(sLattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  // Evaluation
  if (iT % statIter == 0 || converged) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));

    const T L = geom.extentVTI[cfg.flowDirection];

    // probe planes at 10% and 90% (in the structure)
    const T posIn  = T{0.1};
    const T posOut = T{0.9};
    const T Lprobe = (posOut - posIn) * L;

    Vector<T,3> originIn, originOut, extentTK(T{0},T{0},T{0});
    extentTK[cfg.flowDirection] = T{1};

    for (int d = 0; d < 3; ++d) {
      if (d == cfg.flowDirection) {
        originIn[d]  = geom.originVTI[d] + posIn  * geom.extentVTI[d];
        originOut[d] = geom.originVTI[d] + posOut * geom.extentVTI[d];
      } else {
        originIn[d]  = geom.originVTI[d] + T{0.5} * geom.extentVTI[d];
        originOut[d] = geom.originVTI[d] + T{0.5} * geom.extentVTI[d];
      }
    }

    std::vector<int> materials = {1,5};

    // Flux Q through inlet plane
    SuperPlaneIntegralFluxVelocity3D<T> flux(
        sLattice, converter, superGeometry, originIn, extentTK, materials,
        BlockDataReductionMode::Discrete);

    T fluxOut[5]; int dummy[1] {};
    flux(fluxOut, dummy);
    const T Q = fluxOut[0];

    // Average pressures on the two planes
    SuperPlaneIntegralF3D<T> pIn(
        pressure, superGeometry, originIn, extentTK, materials,
        BlockDataReductionMode::Discrete);

    SuperPlaneIntegralF3D<T> pOut(
        pressure, superGeometry, originOut, extentTK, materials,
        BlockDataReductionMode::Discrete);

    T pInArr[2], pOutArr[2];
    pIn(pInArr, dummy);
    pOut(pOutArr, dummy);

    const T p_in   = pInArr[0]  / pInArr[1];
    const T p_out  = pOutArr[0] / pOutArr[1];
    const T deltaP = p_in - p_out;

    // Compute cross-section area orthogonal to flow
    T A = T{1};
    for (int d = 0; d < 3; ++d) {
      if (d != cfg.flowDirection) {
        A *= geom.extentVTI[d];
      }
    }

    results.permeability = (Q * Lprobe * cfg.dynamicViscosity()) / (A * deltaP);

    clout << "----------------------------------------" << std::endl;
    clout << "Flux Q             = " << Q << " m^3/s" << std::endl;
    clout << "deltaP             = " << deltaP << " Pa" << std::endl;
    clout << "Permeability k     = " << results.permeability << " m^2" << std::endl;
    clout << "----------------------------------------" << std::endl;

    // Print the result into permeability.dat
    std::ofstream outfile("permeability.dat");
    outfile << std::scientific << std::setprecision(15) << results.permeability;
    outfile.close();
  }

  // VTK output
  if (converged || iT % vtkIter == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeExternalScalarField3D<T, DESCRIPTOR, descriptors::K> wallK(sLattice);
    wallK.getName() = "wallPermeability";

    vtmWriter.addFunctor(wallK);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.write(iT);
  }
}

int main(int argc, char* argv[])
{
  // === 1) Initialization & command line ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  SimulationConfig cfg;
  std::vector<std::string> cmdInput;

  if (argc > 1) {
    cmdInput.assign(argv + 1, argv + argc);
    cfg.vtiFile = cmdInput[0];
  }
  if (argc > 2)  cfg.arrayName          = cmdInput[1];
  if (argc > 3)  cfg.scalingFactor      = std::stod(cmdInput[2]);
  if (argc > 4)  cfg.uPhys              = std::stod(cmdInput[3]); // physical outlet velocity (m/s)
  if (argc > 5)  cfg.resolution         = std::stoi(cmdInput[4]);
  if (argc > 6)  cfg.InletPressure       = std::stod(cmdInput[5]); // Pa/m
  if (argc > 7)  cfg.tau                = std::stod(cmdInput[6]);
  if (argc > 8)  cfg.wallPermeability   = std::stod(cmdInput[7]);  // m^2, direct target K
  if (argc > 9) cfg.kinematicViscosity = std::stod(cmdInput[8]);   // m^2/s
  if (argc > 10) cfg.fluidDensity       = std::stod(cmdInput[9]);  // kg/m^3
  if (argc > 11) cfg.criterium          = std::stod(cmdInput[10]); 
  if (argc > 12) {
    cfg.flowDirection = std::stoi(cmdInput[11]);
    if (cfg.flowDirection < 0 || cfg.flowDirection > 2) {
      clout << "Error: Flow direction must be 0 (X), 1 (Y), or 2 (Z)" << std::endl;
      return 1;
    }
  } else {
    clout << "Usage: " << argv[0]
          << " <filename> <arrayname> <scaling-factor> <uPhys(m/s)> "
          << "<resolution> <inlet pressure> <tau> <wallPermeability(m^2)> "
          << "<kinematicViscosity(m^2/s)> <fluidDensity(kg/m^3)> <relative tolerance> <flowDirection> "
          << "[uniformGuoZhao(0/1), default 0]\n";
    return 1;
  }
  if (argc > 13) {
    cfg.uniformGuoZhao = (std::stoi(cmdInput[12]) != 0);
  }

  if (cfg.uPhys <= T{0}) {
    clout << "Error: uPhys must be > 0 (m/s)" << std::endl;
    return 1;
  }
  singleton::directories().setOutputDir("./tmp/");


  // Set the directional velocity vector (physical magnitude along flow dir)
  cfg.uVector = {T{0}, T{0}, T{0}};
  cfg.uVector[cfg.flowDirection] = cfg.uPhys;

  clout << "========== Input Parameters ==========\n"
        << "VTI file:         " << cfg.vtiFile << "\n"
        << "Array name:       " << cfg.arrayName << "\n"
        << "Scaling factor:   " << cfg.scalingFactor << "\n"
        << "uPhys (outlet):   " << cfg.uPhys << " m/s\n"
        << "Resolution:       " << cfg.resolution << "\n"
        << "Inlet Pressure:    " << cfg.InletPressure << " Pa/m\n"
        << "Tau:              " << cfg.tau << "\n"
        << "wallPermeability: " << cfg.wallPermeability << " m^2\n"
        << "wallEpsilon:      " << cfg.wallEpsilon << "\n"
        << "Kinematic visc:   " << cfg.kinematicViscosity << " m^2/s\n"
        << "Fluid density:    " << cfg.fluidDensity << " kg/m^3\n"
        << "Tolerance:        " << cfg.criterium << "\n"
        << "Flow direction:   " << (cfg.flowDirection == 0 ? "X" : cfg.flowDirection == 1 ? "Y" : "Z") << "\n"
        << "uniformGuoZhao:   " << (cfg.uniformGuoZhao ? "TRUE (diagnostic mode)" : "false (split dynamics, production)")
        << "\n======================================\n";

  GeometryData      geom;
  RuntimeContext    ctx;
  SimulationResults results;

  // === 2) VTI -> indicator and extents ===
  auto wood = generateIndicatorFromVTI(cfg, cfg.flowDirection, geom, ctx);

  clout << "Domain x-length: " << geom.extent[0] << " m\n";
  clout << "Domain y-length: " << geom.extent[1] << " m\n";
  clout << "Domain z-length: " << geom.extent[2] << " m\n";

  // Characteristic length and characteristic velocity (use uPhys directly)
  const T L = geom.extent[cfg.flowDirection];

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int{cfg.resolution},
      (T)cfg.tau,
      (T)L,                    // Char. length in flow direction
      (T)cfg.uPhys,            // **char physical velocity** (m/s)
      (T)cfg.kinematicViscosity,
      (T)cfg.fluidDensity
  );

  // Simulation horizon a multiple of flow-through time
  const T tFlow    = L / cfg.uPhys;
  const T maxPhysT = T{2.0} * tFlow;  // 2x flow-through time is usually plenty

  clout << "Timeframe to be simulated (≈2×L/u): " << maxPhysT << " s\n";

  converter.print();
  converter.write("permeability");

  clout << "Lattice dx: " << converter.getPhysDeltaX() << " m\n";

  // Wall-permeability report (GuoZhaoBGKdynamics: K is set directly, no floor)
  {
    const T h_SP = converter.getPhysDeltaX();
    clout << "Target wall K:       " << cfg.wallPermeability << " m^2\n";
    clout << "Lattice K (K/h^2):   " << cfg.wallPermeability / (h_SP * h_SP) << "\n";
  }

  // === 3) Decomposition & SuperGeometry ===
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif

  IndicatorLayer3D<T> layer(*wood, converter.getPhysDeltaX());
  CuboidDecomposition3D<T> cuboidDecomposition(layer, converter.getPhysDeltaX(), noOfCuboids);

  // Periodicity: non-flow axes are periodic
  const bool periodicityX = (cfg.flowDirection != 0);
  const bool periodicityY = (cfg.flowDirection != 1);
  const bool periodicityZ = (cfg.flowDirection != 2);
  cuboidDecomposition.setPeriodicity({periodicityX, periodicityY, periodicityZ});

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(cfg, converter, *wood, superGeometry);

  // === 4) Lattice and main loop ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(cfg, converter, *wood, sLattice, superGeometry);

  const int iTmax      = converter.getLatticeTime(maxPhysT);
  const int iTmaxStart = converter.getLatticeTime(maxPhysT * T{0.005});   // ramp over 0.5% of horizon
  bool converged = false;
  std::size_t iT_test = 0;  // Initialize counter
  const std::size_t window = 100;  // Define window size (adjust as needed)
  
  T perm_old = 0.0;

  clout << "Starting simulation...\n";
  clout << "MaxIT: " << iTmax << std::endl;

  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();
  // The main loop for the simulations
  for (std::size_t iT = 0; iT < static_cast<std::size_t>(iTmax); ++iT) {
    setBoundaryValues(sLattice, converter, cfg, geom, static_cast<int>(iT), iTmaxStart, superGeometry);
    sLattice.collideAndStream();
    getResults(sLattice, converter, cfg, geom, results, static_cast<int>(iT), maxPhysT,
               timer, converged, iTmaxStart, superGeometry);

	if ( iT > static_cast<std::size_t>(iTmaxStart) && iT_test >= window ){
		if (fabs(perm_old - results.permeability)/fabs(results.permeability) < cfg.criterium) {
			converged = true;
		getResults(sLattice, converter, cfg, geom, results, static_cast<int>(iT), maxPhysT,
					timer, converged, iTmaxStart, superGeometry);
		clout << "Simulation converged at iT = " << iT << std::endl;
		break;
		}
	   perm_old = results.permeability;
	   iT_test = 0;
	}
	iT_test++;
  }

  timer.stop();
  timer.printSummary();

  return 0;
}
