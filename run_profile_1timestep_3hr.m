%% RUN_PROFILE_1TIMESTEP
% Profiling wrapper: runs exactly 1 timestep of the base case (same config
% as run_input_full2D_pressure2.m, but using item 2's corrected mesh files
% -- solid_leukocyte_P600.mat / solid_endothelium_P300.mat, the version
% pushed 9/17 5:15pm and now the standard baseline) under MATLAB's built-in
% profiler, to find the true computational hot path before parallelizing.
%
% Outputs (written to the working directory):
%   profile_report_1timestep_3hr/   HTML profiler report (profsave)
%   profile_data_1timestep_3hr.mat  raw profile struct, for later inspection
%   out_profile_1timestep_3hr.mat   solver output (out struct)

clc;
clear all;

%% 1. Paths
softlubeDir = '';
if isempty(softlubeDir)
    softlubeDir = fileparts(mfilename('fullpath'));
end
if ~isempty(softlubeDir)
    addpath(softlubeDir);
end

%% 2. Time Controls -- single timestep only
dt = 3e-4;
nSteps = 1;
tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 1024;

%% 3. Global Domain and Pressure Boundary Conditions
zMin = -6e-6;
zMax = 10e-6;
rOuter = 4e-6;

pAtZMin = 0;
pAtZMax = 0;
axisPressureBC = 'dPdr=0';

%% 4. Fluid Mesh Resolution
gap1DWindow = [-0.2e-6, 4.2e-6];
finePressureWindow = gap1DWindow;
coarseDz = 1.0e-6;
fineDz = 0.05e-6;

NrExterior2D = 20;
NrFluid2D = NrExterior2D;
radialFineWindow = [2e-6, 4e-6];
radialFineWeight = 4.0;

%% 5. Plotting and Output -- headless, no figures
plotNative2DPressure = false;
plotNative2DStress = false;
plotNative2DVelocity = false;
plotFirstStepHybridMesh = false;
plotGlobalDomainSchematic = false;
makeLegacyPlots = false;

saveOutput = false;
outputFile = fullfile(softlubeDir, 'simulation_output_profile_1timestep.mat');

closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case (item 2 corrected mesh files)
cfg = struct();

cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endothelium_P300.mat');
cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leukocyte_P600.mat');
cfg.geometry.usePrestressedLeukocyteIC = true;
cfg.geometry.useLeukocyteCenterline = true;

cfg.fluid.zMin = zMin;
cfg.fluid.zMax = zMax;
cfg.fluid.mu = 1.2e-3;
cfg.fluid.slipE = 0e-9;
cfg.fluid.slipL = 0e-9;
cfg.fluid.NrFluid2D = NrFluid2D;
cfg.fluid.Nr = NrFluid2D;
cfg.fluid.NrExterior2D = NrExterior2D;
cfg.fluid.bodyFittedRadialFineWindow = radialFineWindow;
cfg.fluid.bodyFittedRadialFineWeight = radialFineWeight;
cfg.fluid.innerBoundary = 'deformed_leukocyte';
cfg.fluid.useSlopeAwareWallKinematics = true;

cfg.bc.pIn = pAtZMin;
cfg.bc.pOut = pAtZMax;
cfg.bc.UwL = 0;
cfg.bc.UwE = 0;

cfg.solid.noLeukocyte = false;
cfg.solid.leukocyte.enabled = true;
cfg.solid.leukocyte.deformable = true;
cfg.solid.leukocyte.EL = 200;
cfg.solid.leukocyte.nuL = 0.46;
cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
cfg.solid.leukocyte.etaL = 1;
cfg.solid.leukocyte.supportL = 'axis';

cfg.solid.endothelium.Ee = 500;
cfg.solid.endothelium.nuE = 0.46;
cfg.solid.endothelium.useViscoelasticEndothelium = true;
cfg.solid.endothelium.etaE = 1;

cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;
cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
cfg.interface.warnPrestressLoadMismatch = false;
cfg.numerics.dt = dt;
cfg.numerics.tEnd = tEnd;
cfg.numerics.minGap = 0.001e-6;
cfg.numerics.maxNewtonMono = 100;
cfg.numerics.tolNewtonMono = 1e-6;
cfg.numerics.monoSolidAbsTol = 1e-7;
cfg.numerics.monoFluidAbsTol = 1e-8;
cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
cfg.numerics.enablePressureJumpRetry = false;
cfg.numerics.maxPressureJumpAbs = inf;
cfg.numerics.maxPressureJump2DAbs = inf;
cfg.numerics.useObjectiveKelvinVoigt = true;
cfg.numerics.useDirectFluidSolve = true;
cfg.numerics.useFsolveMono = true;
cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
cfg.numerics.dtMin = dtMin;
cfg.numerics.dtRetryFactor = 0.5;
cfg.numerics.dtGrowFactor = 1.25;
cfg.numerics.maxTimeStepRetries = 18;
cfg.numerics.fsolveDisplay = 'off';
cfg.numerics.diagnosticsEnabled = true;

cfg.output.saveOutput = saveOutput;
cfg.output.outputFile = outputFile;
cfg.output.makePlots = makeLegacyPlots;
cfg.output.plotFinalGlobalPressureContour = false;
cfg.output.plotGlobal2DPressureTractionComparison = false;
cfg.output.printEvery = printEvery;
cfg.output.storeFull2DFluidHist = false;
cfg.output.store2DFluidArrays = true;
cfg.output.store2DPhysicalGridHist = true;
cfg.output.store2DSpeedHist = true;

cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = closeFiguresAtStart;
cfg.parOverrides = struct();

%% 7. Global Axial Mesh and Exterior Boundary Embedding
cfg.global1D.zDomain = [zMin, zMax];
cfg.global1D.rDomain = [0, rOuter];
cfg.global1D.axisPressureBC = axisPressureBC;
cfg.global1D.outerPressureBC = 0;

cfg.parOverrides.useGlobal1DPressure = true;
cfg.parOverrides.global1DOuterRadius = rOuter;
cfg.parOverrides.global1DAxisRadius = 1e-9;
cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
cfg.parOverrides.global1DFineWindow = finePressureWindow;
cfg.parOverrides.global1DCoarseDz = coarseDz;
cfg.parOverrides.global1DFineDz = fineDz;
cfg.parOverrides.global1DPlotNr = 161;

%% 8. Hybrid Gap-1D / Exterior-2D Fluid Module
cfg.fluid.useFull2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = true;
cfg.fluid.hybridGapZ = gap1DWindow;
cfg.fluid.hybridExteriorMinCells = 4;
cfg.fluid.hybridExteriorFailMode = 'warn';
cfg.fluid.hybridUseExterior2DPressureVector = true;
cfg.fluid.useBodyFittedMACFluid = true;
cfg.fluid.useBodyFittedMACTractionCorrection = true;
cfg.fluid.useFeedbackTractionCorrection = false;
cfg.fluid.useAitkenTractionCorrectionRelax = false;
cfg.fluid.debugRadialAxialTractionCorrection = false;
cfg.fluid.maxBodyFittedTractionCorrections = 100;
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.5;
cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.fluid.resolveFluidAfterTractionCorrection = true;
cfg.fluid.useBodyFittedMACTractionInSolid = true;

cfg.parOverrides.useFull2DFluid = cfg.fluid.useFull2DFluid;
cfg.parOverrides.useHybridGap1DExterior2DFluid = cfg.fluid.useHybridGap1DExterior2DFluid;
cfg.parOverrides.hybridGapZ = cfg.fluid.hybridGapZ;
cfg.parOverrides.NrExterior2D = NrExterior2D;
cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
cfg.parOverrides.hybridExteriorMinCells = cfg.fluid.hybridExteriorMinCells;
cfg.parOverrides.hybridExteriorFailMode = cfg.fluid.hybridExteriorFailMode;
cfg.parOverrides.hybridUseExterior2DPressureVector = cfg.fluid.hybridUseExterior2DPressureVector;
cfg.parOverrides.useBodyFittedMACFluid = cfg.fluid.useBodyFittedMACFluid;
cfg.parOverrides.useBodyFittedMACTractionCorrection = cfg.fluid.useBodyFittedMACTractionCorrection;
cfg.parOverrides.useFeedbackTractionCorrection = cfg.fluid.useFeedbackTractionCorrection;
cfg.parOverrides.useAitkenTractionCorrectionRelax = cfg.fluid.useAitkenTractionCorrectionRelax;
cfg.parOverrides.maxBodyFittedTractionCorrections = cfg.fluid.maxBodyFittedTractionCorrections;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = cfg.fluid.bodyFittedTractionCorrectionRelax;
cfg.parOverrides.bodyFittedTractionCorrectionFailMode = cfg.fluid.bodyFittedTractionCorrectionFailMode;
cfg.parOverrides.resolveFluidAfterTractionCorrection = cfg.fluid.resolveFluidAfterTractionCorrection;
cfg.parOverrides.useBodyFittedMACTractionInSolid = cfg.fluid.useBodyFittedMACTractionInSolid;
cfg.parOverrides.NrFluid2D = NrFluid2D;
cfg.parOverrides.Nr = NrFluid2D;
cfg.parOverrides.useGlobal2DPressureTraction = false;

%% 9. Run under the profiler
fprintf('\nProfiling 1-timestep base case (item 2 corrected mesh files)\n');
fprintf('   dt                  = %.6e s\n', dt);
fprintf('   tEnd                = %.6e s\n', tEnd);
fprintf('   requested steps     = %.0f\n\n', nSteps);

profile('clear');
profile('on');

tStart = tic;
out = softlube_run_case_global_coupled(cfg);
wallClockSeconds = toc(tStart);

profile('off');

fprintf('\n1-timestep profiling run finished\n');
fprintf('   stopStep       = %d\n', out.stopStep);
fprintf('   final t        = %.6e s\n', out.t(end));
fprintf('   wall-clock     = %.3f s\n', wallClockSeconds);

pInfo = profile('info');
save('profile_data_1timestep_3hr.mat', 'pInfo', 'wallClockSeconds');
profsave(pInfo, 'profile_report_1timestep_3hr');

out.cfg = cfg;
save('out_profile_1timestep_3hr.mat', 'out', 'wallClockSeconds');

fprintf('\nSaved: profile_data_1timestep_3hr.mat, profile_report_1timestep_3hr/, out_profile_1timestep_3hr.mat\n');
