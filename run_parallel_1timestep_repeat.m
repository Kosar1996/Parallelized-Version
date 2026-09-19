%% RUN_PARALLEL_1TIMESTEP_REPEAT
% Reproducibility diagnostic: identical to run_parallel_1timestep.m, but
% saves to out_parallel_1timestep_repeat_N<P>.mat instead of overwriting the
% original result. Run with the SAME NUM_PROCS as an earlier run (typically
% N=1) to check whether the solver's final converged state is even
% reproducible for a fixed processor count -- the first-round results showed
% max|P2D| varying wildly across different N (426 Pa at N=1, up to 2533 Pa
% at N=8), and this test isolates whether that's really tied to processor
% count or is a more general non-determinism.
%
% Outputs (written to the working directory, suffixed by processor count):
%   out_parallel_1timestep_repeat_N<P>.mat   solver output (out struct) + timing

clc;
clear all;

%% 0. Processor count and parallel pool
numProcsEnv = str2double(getenv('NUM_PROCS'));
if ~isfinite(numProcsEnv) || numProcsEnv < 1
    error(['NUM_PROCS environment variable must be set to a positive integer ', ...
        'before running this script (the submit script sets it).']);
end
numProcs = round(numProcsEnv);

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool('local', numProcs);

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
outputFile = fullfile(softlubeDir, sprintf('simulation_output_parallel_1timestep_repeat_N%d.mat', numProcs));

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

%% 9. Run
fprintf('\nParallel 1-timestep base case, NUM_PROCS = %d\n', numProcs);
fprintf('   dt                  = %.6e s\n', dt);
fprintf('   tEnd                = %.6e s\n', tEnd);
fprintf('   requested steps     = %.0f\n\n', nSteps);

tStart = tic;
out = softlube_run_case_global_coupled(cfg);
wallClockSeconds = toc(tStart);

fprintf('\n1-timestep parallel run finished (NUM_PROCS = %d)\n', numProcs);
fprintf('   stopStep       = %d\n', out.stopStep);
fprintf('   final t        = %.6e s\n', out.t(end));
fprintf('   wall-clock     = %.3f s\n', wallClockSeconds);

out.cfg = cfg;
out.numProcs = numProcs;
outFileName = sprintf('out_parallel_1timestep_repeat_N%d.mat', numProcs);
save(outFileName, 'out', 'wallClockSeconds', 'numProcs');

fprintf('\nSaved: %s\n', outFileName);

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
