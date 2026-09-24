%% RUN_PARALLEL_150_WIDE6
% 150-timestep parallel benchmark on the new codebase package: full 2D
% fluid solver (useFull2DFluid=true, useHybridGap1DExterior2DFluid=false),
% leukocyte prestress P300 (was P600), outer domain radius 6um (was 4um).
% Base config matches run_case_7_300_wide6.m exactly; nSteps overridden to
% 150 per the benchmark spec (baseline: <6 wall-clock hours at N=1).
%
% Same parallel-pool pattern as the earlier 390-step scripts (parpool
% opened for every processor count including N=1, for a fair comparison,
% main-thread pinned to 1 computation thread for run-to-run
% reproducibility). Uses the parfor-restructured assembly functions
% (assemble_finite_def_axisym.m and the other 3), already verified
% bit-identical across processor counts on both the pre-fix and fixed-
% physics 390-timestep runs -- the underlying physics in these 4
% functions is unchanged in this new codebase (confirmed byte-identical
% to the known serial baseline before reapplying parfor).
%
% Outputs (written to the working directory, suffixed by processor count):
%   out_parallel_150_wide6_N<P>.mat   solver output (out struct) + timing

clc;
clear all;

%% 0. Processor count and parallel pool
numProcsEnv = str2double(getenv('NUM_PROCS'));
if ~isfinite(numProcsEnv) || numProcsEnv < 1
    error(['NUM_PROCS environment variable must be set to a positive integer ', ...
        'before running this script (the submit script sets it).']);
end
numProcs = round(numProcsEnv);

try
    maxNumCompThreads(1);
catch ME
    warning('maxNumCompThreads not available/settable: %s', ME.message);
end

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

%% 2. CUSTOMIZABLE PARAMETERS (matches run_case_7_300_wide6.m)
dt = 1e-5;
nSteps = 150;

leukocyteDeformable = true;
leukocyteUseViscoelastic = true;
leukocyteEtaL = 1;

endotheliumUseViscoelastic = true;
endotheliumEtaE = 1;

tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 1024;

%% 3. Global Domain and Pressure Boundary Conditions
zMin = -6e-6;
zMax = 10e-6;
rOuter = 6e-6;

pAtZMin = 0;
pAtZMax = 0;
axisPressureBC = 'dPdr=0';

%% 4. Fluid Mesh Resolution
gap1DWindow = [-0.5e-6, 4.5e-6];
finePressureWindow = gap1DWindow;
coarseDz = 1.0e-6;
fineDz = 0.05e-6;

NrExterior2D = 20;
NrFluid2D = NrExterior2D;
radialFineWindow = [2e-6, rOuter];
radialFineWeight = 4.0;

%% 5. Plotting and Output -- headless, no figures
saveOutput = false;
outputFile = fullfile(softlubeDir, sprintf('simulation_output_150_wide6_N%d.mat', numProcs));
closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case
cfg = struct();

cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endothelium_P300.mat');
cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leukocyte_P300.mat');
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
cfg.solid.leukocyte.deformable = leukocyteDeformable;
cfg.solid.leukocyte.useViscoelasticLeukocyte = leukocyteUseViscoelastic;
cfg.solid.leukocyte.etaL = leukocyteEtaL;
cfg.solid.leukocyte.supportL = 'axis';

cfg.solid.endothelium.useViscoelasticEndothelium = endotheliumUseViscoelastic;
cfg.solid.endothelium.etaE = endotheliumEtaE;

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
cfg.output.checkpointFile = outputFile;
cfg.output.checkpointEvery = 1;
cfg.output.makePlots = false;
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

%% 8. Full 2D Fluid Module (hybrid disabled -- full 2D solver)
cfg.fluid.useFull2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = false;
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

cfg.fluid.hybridTransitionBuffer = 0.3e-6;
cfg.fluid.useSmoothHybridBlending = true;

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
cfg.parOverrides.useGlobal2DPressureTraction = true;

cfg.parOverrides.hybridTransitionBuffer = cfg.fluid.hybridTransitionBuffer;
cfg.parOverrides.useSmoothHybridBlending = cfg.fluid.useSmoothHybridBlending;

%% 9. Run
fprintf('\nParallel 150-timestep wide6/P300 benchmark, NUM_PROCS = %d\n', numProcs);
fprintf('   dt   = %.6e s\n', dt);
fprintf('   tEnd = %.6e s\n\n', tEnd);

tStart = tic;
out = softlube_run_case_global_coupled(cfg);
wallClockSeconds = toc(tStart);

fprintf('\n150-timestep benchmark finished (NUM_PROCS = %d)\n', numProcs);
fprintf('   stopStep       = %d\n', out.stopStep);
fprintf('   final t        = %.6e s\n', out.t(end));
fprintf('   wall-clock     = %.3f s (%.2f hr)\n', wallClockSeconds, wallClockSeconds/3600);

out.cfg = cfg;
out.numProcs = numProcs;
outFileName = sprintf('out_parallel_150_wide6_N%d.mat', numProcs);
save(outFileName, 'out', 'wallClockSeconds', 'numProcs', '-v7.3');

fprintf('\nSaved: %s\n', outFileName);

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
