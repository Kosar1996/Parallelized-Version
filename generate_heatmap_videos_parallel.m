function generate_heatmap_videos_parallel(matFile, label, outDir, fps)
% GENERATE_HEATMAP_VIDEOS_PARALLEL Multi-core rendering pipeline for 2D heatmaps.
% Uses repmat across 3D array slice dimensions to eliminate zero-padding domain
% truncations, fixing both fluid voids and stress plot overlays.

if nargin < 4 || isempty(fps)
    fps = 5;
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Suppress non-fatal headless graphics acceleration warning on HPC client
warning('off', 'MATLAB:graphics:noGraphicsAcceleration');

% Force software OpenGL rendering for headless cluster stability
try
    opengl('save', 'software');
catch
end

% Ensure Parallel Pool is Active
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting MATLAB parallel pool...\n');
    parpool(); 
end

% Broadcast warning suppressions to all active worker processes
pctRunOnAll warning('off', 'MATLAB:graphics:noGraphicsAcceleration');
pctRunOnAll warning('off', 'MATLAB:structOnObject');

S = load(matFile);
if isfield(S, 'out')
    out = S.out;
elseif isfield(S, 'outRestart')
    out = S.outRestart;
else
    error('generate_heatmap_videos_parallel:noOut', 'No ''out'' or ''outRestart'' found in %s', matFile);
end

nSteps = out.stopStep;
if nSteps < 1
    error('generate_heatmap_videos_parallel:noSteps', '%s has stopStep=%d', matFile, nSteps);
end

fields = {'pressure', @plot_select_native2d_pressure; ...
          'stress',   @plot_select_native2d_stress; ...
          'velocity', @plot_select_native2d_velocity};

% Create temporary directory for worker-rendered image frames
tempFrameDir = fullfile(outDir, ['temp_frames_' label]);
if ~exist(tempFrameDir, 'dir')
    mkdir(tempFrameDir);
end

for f = 1:size(fields, 1)
    fieldName = fields{f,1};
    plotFn    = fields{f,2};
    fprintf('\n=== Parallel Rendering [%s]: %s ===\n', label, fieldName);

    fieldTempDir = fullfile(tempFrameDir, fieldName);
    if ~exist(fieldTempDir, 'dir')
        mkdir(fieldTempDir);
    end

    % ---------------------------------------------------------------------
    % Step A: Slice History Arrays & Build Base Context
    % ---------------------------------------------------------------------
    stateHist_sliced  = out.stateHist(1:nSteps);
    fluidHist_sliced  = out.fluidHist(1:nSteps);
    PHist_sliced      = out.PHist(:,:,1:nSteps);
    urCHist_sliced    = out.urCHist(:,:,1:nSteps);
    uzCHist_sliced    = out.uzCHist(:,:,1:nSteps);
    speedCHist_sliced = out.speedCHist(:,:,1:nSteps);

    % Build complete base structure
    outBase = struct();
    outBase.meshE      = out.meshE;      % Endothelial solid mesh
    outBase.meshL      = out.meshL;      % Leukocyte solid mesh
    outBase.par        = out.par;        % Model parameters
    outBase.t          = out.t;          % Full time vector
    outBase.stopStep   = nSteps;         % Total step count
    
    % Pass complete fluidHist cell array required for boundary masking
    outBase.fluidHist  = out.fluidHist(1:nSteps);

    % Copy interface and grid definitions required for domain masking
    if isfield(out, 'meshF'),      outBase.meshF      = out.meshF; end
    if isfield(out, 'interfaceE'), outBase.interfaceE = out.interfaceE; end
    if isfield(out, 'interfaceL'), outBase.interfaceL = out.interfaceL; end
    if isfield(out, 'z'),          outBase.z          = out.z; end
    if isfield(out, 'zGrid'),      outBase.zGrid      = out.zGrid; end
    if isfield(out, 'dz'),         outBase.dz         = out.dz; end
    if isfield(out, 'dtHist'),     outBase.dtHist     = out.dtHist; end

    % Render t=0 reference frame sequentially if available
    hasT0Frame = isfield(out, 't0State') && isfield(out, 't0Fluid') && ~isempty(out.t0Fluid);
    if hasT0Frame
        try
            out0 = outBase;
            out0.fluidHist = {out.t0Fluid};
            out0.stateHist = {out.t0State};
            out0.state     = out.t0State;
            out0.fluid     = out.t0Fluid;
            if isfield(out.t0Fluid, 'meshF')
                out0.meshF = out.t0Fluid.meshF;
            end
            out0.stopStep  = 1;
            out0.t         = 0;
            out0.PHist     = out.t0Fluid.P;
            out0.native2D  = struct();
            out0.native2D.P     = out.t0Fluid.P;
            out0.native2D.R     = out.t0Fluid.meshF.Rp;
            out0.native2D.Z     = out.t0Fluid.meshF.Zp;
            out0.native2D.ur    = out.t0Fluid.urC;
            out0.native2D.uz    = out.t0Fluid.uzC;
            out0.native2D.speed = sqrt(out.t0Fluid.urC.^2 + out.t0Fluid.uzC.^2);
            out0.useHybridGap1DExterior2DFluid = true;
            out0.useFull2DFluid = true;

            fig0 = figure('Visible', 'off');
            plotFn(out0, 1);
            
            frame0 = getframe(fig0);
            imwrite(frame0.cdata, fullfile(fieldTempDir, 'frame_000000.png'));
            close(fig0);
            delete(fig0);
        catch ME0
            fprintf('  t=0 frame failed: %s\n', ME0.message);
        end
    end

    % ---------------------------------------------------------------------
    % Step B: Parallel Frame Generation (parfor)
    % ---------------------------------------------------------------------
    parfor k = 1:nSteps
        fig = [];
        try
            outK = outBase;
            
            % Assign step k into stateHist and fluidHist cell arrays
            outK.stateHist = cell(1, nSteps);
            outK.stateHist{k} = stateHist_sliced{k};
            outK.state = stateHist_sliced{k};

            outK.fluid = fluidHist_sliced{k};

            % Carry displacement fields for solid stress recovery
            if isfield(stateHist_sliced{k}, 'uE')
                outK.state.uE = stateHist_sliced{k}.uE;
            end
            if isfield(stateHist_sliced{k}, 'uL')
                outK.state.uL = stateHist_sliced{k}.uL;
            end

            % Ensure fluid mesh coordinates carry current step Rp/Zp definitions
            if isfield(outK.fluid, 'meshF')
                outK.meshF = outK.fluid.meshF;
            end

            % Force hybrid field recognition
            outK.useHybridGap1DExterior2DFluid = true;
            outK.useFull2DFluid = true;

            % Extract step k 2D matrices
            Pk   = PHist_sliced(:,:,k);
            Rk   = out.RPHist(:,:,k);
            Zk   = out.ZPHist(:,:,k);
            urk  = urCHist_sliced(:,:,k);
            uzk  = uzCHist_sliced(:,:,k);
            spdk = speedCHist_sliced(:,:,k);

            % REPMAT step k 2D slice across 3rd dimension to eliminate zeros
            outK.PHist   = repmat(Pk,  [1, 1, nSteps]);
            outK.RPHist  = repmat(Rk,  [1, 1, nSteps]);
            outK.ZPHist  = repmat(Zk,  [1, 1, nSteps]);
            outK.urCHist = repmat(urk, [1, 1, nSteps]);
            outK.uzCHist = repmat(uzk, [1, 1, nSteps]);

            % Populate native2D structure for fluid stress overlay
            outK.native2D = struct();
            outK.native2D.P     = Pk;
            outK.native2D.R     = Rk;
            outK.native2D.Z     = Zk;
            outK.native2D.ur    = urk;
            outK.native2D.uz    = uzk;
            outK.native2D.speed = spdk;

            fig = figure('Visible', 'off');
            
            % Render plot at step k
            plotFn(outK, k);
            
            frame = getframe(fig);
            imgFile = fullfile(fieldTempDir, sprintf('frame_%06d.png', k));
            imwrite(frame.cdata, imgFile);
            
            close(fig);
            delete(fig);
        catch MEk
            fprintf('  Worker frame %d failed: %s\n', k, MEk.message);
            if ~isempty(fig) && ishandle(fig)
                close(fig);
                delete(fig);
            end
        end
    end

    % ---------------------------------------------------------------------
    % Step C: Sequential Video Assembly
    % ---------------------------------------------------------------------
    imgFiles = dir(fullfile(fieldTempDir, 'frame_*.png'));
    if isempty(imgFiles)
        fprintf('  [Warning] No frames rendered for %s; skipping video output.\n', fieldName);
        continue;
    end

    try
        outFile = fullfile(outDir, [label '_' fieldName '.mp4']);
        v = VideoWriter(outFile, 'MPEG-4');
    catch
        outFile = fullfile(outDir, [label '_' fieldName '.avi']);
        v = VideoWriter(outFile, 'Motion JPEG AVI');
    end
    v.FrameRate = fps;

    open(v);
    for idx = 1:numel(imgFiles)
        img = imread(fullfile(fieldTempDir, imgFiles(idx).name));
        writeVideo(v, img);
    end
    close(v);
    
    fprintf('  Completed: %d frames compiled -> %s\n', numel(imgFiles), outFile);
end

% Clean up temporary PNG working directory
if exist(tempFrameDir, 'dir')
    rmdir(tempFrameDir, 's');
end

fprintf('\nParallel video generation completed successfully for: %s\n', label);
end
