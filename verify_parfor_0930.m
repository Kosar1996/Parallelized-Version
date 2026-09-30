%% VERIFY_PARFOR_0930
% Unit check for the 0930 parallelization, run under a REAL local parallel
% pool (parfor classification problems only surface with an open pool).
% Compares each of the 4 parallelized assembly functions in this folder
% against the untouched serial versions in ../0930_serial (= the
% Leukocyte_Main_Files-0928_v2 package) on the production endothelium mesh.
% PASS criterion: outputs EXACTLY identical (isequal, difference = 0),
% not merely within a tolerance.

clc;
clear all;
parDir    = fileparts(mfilename('fullpath'));
serialDir = fullfile(parDir, '..', '0930_serial');
if ~exist(fullfile(serialDir, 'assemble_finite_def_axisym.m'), 'file')
    error('Serial folder not found at %s (upload 0930_serial next to 0930_parallel).', serialDir);
end

try, maxNumCompThreads(1); catch, end
p = gcp('nocreate'); if ~isempty(p), delete(p); end
parpool('local', 2);

% IMPORTANT: run from a neutral folder. MATLAB always prefers files in the
% current folder over the path, so staying in parDir would make the
% "serial" reference call silently use the parallel files.
cd(tempdir);
addpath(parDir, '-begin');
S = load(fullfile(parDir, 'solid_endothelium_P300.mat'));
par = S.par;
par.useViscoelasticEndothelium = true;
par.etaE = 1;
par.dt = 2.5e-6;
par.useObjectiveKelvinVoigt = true;
meshRaw = S.meshE;
meshC   = prepare_axisym_mesh_cache(meshRaw);

ndof = size(meshRaw.nodes,1) * 2;
% Physical test state: the prestressed endothelium displacement from the
% production initialization file (what step 1 of the benchmark starts
% from), plus a tiny previous-state offset so the viscous terms are nonzero.
if isfield(S, 'uE_pre')
    u = S.uE_pre(:);
else
    u = zeros(ndof, 1);
end
rng(0);
uOld = u + 1e-13 * randn(ndof, 1);

allPass = true;
cases = {'cached', meshC; 'uncached', meshRaw};

for c = 1:size(cases,1)
    tagc = cases{c,1}; mesh = cases{c,2};
    fprintf('\n===== mesh: %s =====\n', tagc);

    % --- serial reference (from ../0930_serial) ---
    addpath(serialDir, '-begin'); rehash path;
    fprintf('  serial reference from: %s\n', which('assemble_finite_def_axisym'));
    [F1s, K1s] = assemble_finite_def_axisym(mesh, u, par);
    F2s = assemble_finite_def_internal_force_only(mesh, u, par);
    F4s = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par);
    if strcmp(tagc, 'cached')
        [F3s, K3s] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
    end
    rmpath(serialDir); addpath(parDir, '-begin'); rehash path;

    % --- parallel (this folder) ---
    fprintf('  parallel version from: %s\n', which('assemble_finite_def_axisym'));
    [F1p, K1p] = assemble_finite_def_axisym(mesh, u, par);
    F2p = assemble_finite_def_internal_force_only(mesh, u, par);
    F4p = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par);
    if strcmp(tagc, 'cached')
        [F3p, K3p] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
    end

    allPass = report('assemble_finite_def_axisym  Fint', F1s, F1p) && allPass;
    allPass = report('assemble_finite_def_axisym  K',    K1s, K1p) && allPass;
    allPass = report('assemble_finite_def_internal_force_only', F2s, F2p) && allPass;
    allPass = report('assemble_axisym_kelvin_voigt_viscous_force_only', F4s, F4p) && allPass;
    if strcmp(tagc, 'cached')
        allPass = report('assemble_axisym_kelvin_voigt_viscous  Fvisc', F3s, F3p) && allPass;
        allPass = report('assemble_axisym_kelvin_voigt_viscous  Kvisc', K3s, K3p) && allPass;
    else
        fprintf('  (KV tangent, uncached branch: not tested -- that branch calls an undefined\n');
        fprintf('   function in the ORIGINAL serial code too and is never reached in production)\n');
    end
end

fprintf('\n==========================================\n');
if allPass
    fprintf('OVERALL: PASS -- all outputs bit-identical to serial (diff = 0)\n');
else
    fprintf('OVERALL: FAIL -- see above\n');
end
fprintf('==========================================\n');
p = gcp('nocreate'); if ~isempty(p), delete(p); end

function ok = report(name, A, B)
    if issparse(A), d = full(max(abs(A(:) - B(:)))); else, d = max(abs(A(:) - B(:))); end
    if isempty(d), d = 0; end
    ok = isequal(A, B);
    fprintf('  %-52s max|diff| = %.3e   %s\n', name, d, ternary(ok, 'PASS (identical)', 'FAIL'));
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
