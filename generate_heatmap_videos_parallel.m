function generate_heatmap_videos_parallel(matFile, label, outDir, fps)
% GENERATE_HEATMAP_VIDEOS_PARALLEL Multi-core rendering pipeline for 2D heatmaps.
% Enforces fixed pixel figure positions and explicit axis bounds across parfor workers
% to eliminate frame-to-frame canvas resizing and heatmap popping.

if nargin < 4 || isempty(fps)
    fps = 5;
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Suppress non-fatal headless graphics acceleration warning on HPC client
warning('off', 'MATLAB:graphics:noGraphicsAcceleration');

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
    % Step A: Build Slice History Context
    % ---------------------------------------------------------------------
    stateHist_sliced = out.stateHist(1:nSteps);
    fluidHist_sliced = out.fluidHist(1:nSteps);

    if isfield(out, 'PHist') && ~isempty(out.PHist)
        PHist_sliced = out.PHist(:,:,1:nSteps);
    else
        PHist_sliced = zeros(size(fluidHist_sliced{1}.P,1), size(fluidHist_sliced{1}.P,2), nSteps);
        for stIdx = 1:nSteps
            PHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.P;
        end
    end

    if isfield(out, 'RPHist') && ~isempty(out.RPHist)
        RPHist_sliced = out.RPHist(:,:,1:nSteps);
    else
        RPHist_sliced = zeros(size(fluidHist_sliced{1}.meshF.Rp,1), size(fluidHist_sliced{1}.meshF.Rp,2), nSteps);
        for stIdx = 1:nSteps
            RPHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.meshF.Rp;
        end
    end

    if isfield(out, 'ZPHist') && ~isempty(out.ZPHist)
        ZPHist_sliced = out.ZPHist(:,:,1:nSteps);
    else
        ZPHist_sliced = zeros(size(fluidHist_sliced{1}.meshF.Zp,1), size(fluidHist_sliced{1}.meshF.Zp,2), nSteps);
        for stIdx = 1:nSteps
            ZPHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.meshF.Zp;
        end
    end

    if isfield(out, 'urCHist') && ~isempty(out.urCHist)
        urCHist_sliced = out.urCHist(:,:,1:nSteps);
    else
        urCHist_sliced = zeros(size(fluidHist_sliced{1}.urC,1), size(fluidHist_sliced{1}.urC,2), nSteps);
        for stIdx = 1:nSteps
            urCHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.urC;
        end
    end

    if isfield(out, 'uzCHist') && ~isempty(out.uzCHist)
        uzCHist_sliced = out.uzCHist(:,:,1:nSteps);
    else
        uzCHist_sliced = zeros(size(fluidHist_sliced{1}.uzC,1), size(fluidHist_sliced{1}.uzC,2), nSteps);
        for stIdx = 1:nSteps
            uzCHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.uzC;
        end
    end

    outBase = struct();
    outBase.meshE      = out.meshE;      
    outBase.meshL      = out.meshL;      
    outBase.par        = out.par;        
    outBase.t          = out.t;          
    outBase.stopStep   = nSteps;         
    outBase.fluidHist  = fluidHist_sliced;
    outBase.stateHist  = stateHist_sliced;
    outBase.PHist      = PHist_sliced;
    outBase.RPHist     = RPHist_sliced;
    outBase.ZPHist     = ZPHist_sliced;
    outBase.urCHist    = urCHist_sliced;
    outBase.uzCHist    = uzCHist_sliced;

    if isfield(out, 'meshF'),      outBase.meshF      = out.meshF; end
    if isfield(out, 'interfaceE'), outBase.interfaceE = out.interfaceE; end
    if isfield(out, 'interfaceL'), outBase.interfaceL = out.interfaceL; end
    if isfield(out, 'z'),          outBase.z          = out.z; end
    if isfield(out, 'zGrid'),      outBase.zGrid      = out.zGrid; end
    if isfield(out, 'dz'),         outBase.dz         = out.dz; end
    if isfield(out, 'dtHist'),     outBase.dtHist     = out.dtHist; end

    % ---------------------------------------------------------------------
    % Step B: Parallel Frame Generation with Enforced Canvas & Bounds Lock
    % ---------------------------------------------------------------------
    parfor k = 1:nSteps
        fig = [];
        try
            outK = outBase;
            outK.stateHist = cell(1, nSteps);
            outK.stateHist{k} = stateHist_sliced{k};
            outK.state = stateHist_sliced{k};
            outK.fluid = fluidHist_sliced{k};

            if isfield(stateHist_sliced{k}, 'uE')
                outK.state.uE = stateHist_sliced{k}.uE;
            end
            if isfield(stateHist_sliced{k}, 'uL')
                outK.state.uL = stateHist_sliced{k}.uL;
            end

            if isfield(outK.fluid, 'meshF')
                outK.meshF = outK.fluid.meshF;
            end

            outK.useHybridGap1DExterior2DFluid = true;
            outK.useFull2DFluid = true;

            Pk   = PHist_sliced(:,:,k);
            Rk   = RPHist_sliced(:,:,k);
            Zk   = ZPHist_sliced(:,:,k);
            urk  = urCHist_sliced(:,:,k);
            uzk  = uzCHist_sliced(:,:,k);

            outK.PHist   = repmat(Pk,  [1, 1, nSteps]);
            outK.RPHist  = repmat(Rk,  [1, 1, nSteps]);
            outK.ZPHist  = repmat(Zk,  [1, 1, nSteps]);
            outK.urCHist = repmat(urk, [1, 1, nSteps]);
            outK.uzCHist = repmat(uzk, [1, 1, nSteps]);

            outK.native2D = struct();
            outK.native2D.P     = Pk;
            outK.native2D.R     = Rk;
            outK.native2D.Z     = Zk;
            outK.native2D.ur    = urk;
            outK.native2D.uz    = uzk;

            % 1. Instantiate Figure with EXPLICIT Fixed Window Pixel Position
            fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1200 800]);
            plotFn(outK, k);
            
            % 2. Lock Axis Position & Domain Limits
            ax = gca;
            set(ax, 'Units', 'normalized');
            if ~strcmp(fieldName, 'stress')
                % Lock exact plot canvas area (left, bottom, width, height)
                set(ax, 'Position', [0.12 0.12 0.75 0.80], 'ActivePositionProperty', 'position');
                xlim(ax, [0 4]);
                ylim(ax, [-6 10]);
            end

            % 3. Capture Exact Figure Canvas Frame
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
    % Step C: Sequential AVI Assembly + Native FFmpeg MP4 Conversion
    % ---------------------------------------------------------------------
    imgFiles = dir(fullfile(fieldTempDir, 'frame_*.png'));
    if isempty(imgFiles)
        fprintf('  [Warning] No frames rendered for %s; skipping video output.\n', fieldName);
        continue;
    end

    tempAviFile = fullfile(outDir, [label '_' fieldName '_temp.avi']);
    finalMp4File = fullfile(outDir, [label '_' fieldName '.mp4']);

    v = VideoWriter(tempAviFile, 'Motion JPEG AVI');
    v.FrameRate = fps;
    if isprop(v, 'Quality'), v.Quality = 100; end

    open(v);
    for idx = 1:numel(imgFiles)
        img = imread(fullfile(fieldTempDir, imgFiles(idx).name));
        writeVideo(v, img);
    end
    close(v);

    % Native MPEG-4 conversion via FFmpeg using mpeg4 codec
    cmd = sprintf('ffmpeg -y -i "%s" -c:v mpeg4 -q:v 2 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" "%s"', ...
        tempAviFile, finalMp4File);
    [status, cmdOut] = system(cmd);

    if status == 0 && exist(finalMp4File, 'file')
        delete(tempAviFile);
        fprintf('  Successfully compiled smooth MP4: %d frames -> %s\n', numel(imgFiles), finalMp4File);
    else
        actualAviFile = fullfile(outDir, [label '_' fieldName '.avi']);
        movefile(tempAviFile, actualAviFile);
        fprintf('  [Warning] FFmpeg conversion notice (%s). Saved AVI: %s\n', strtrim(cmdOut), actualAviFile);
    end
end

if exist(tempFrameDir, 'dir')
    rmdir(tempFrameDir, 's');
end

fprintf('\nParallel video generation completed successfully for: %s\n', label);
end
