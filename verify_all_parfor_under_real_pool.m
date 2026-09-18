%% VERIFY_ALL_PARFOR_UNDER_REAL_POOL
% Re-verification after the cell-array parfor fix. The earlier verify_*.m
% scripts passed locally but the cluster run still failed, because MATLAB's
% parfor variable-classification check only actually triggers when a real
% parallel pool is open -- without one, parfor silently falls back to a
% serial-like path that never exercises the classification logic. This
% script opens a real local pool first, so it actually tests what the
% cluster run does.
%
% Compares each of the 4 hot functions' NEW (live, cell-array-based parfor)
% version against the pre_parallelization_backup/ (original serial) version,
% under an actual open parpool.

clc;
clear all;
softlubeDir = fileparts(mfilename('fullpath'));
addpath(softlubeDir);
backupDir = fullfile(softlubeDir, 'pre_parallelization_backup');

fprintf('Opening a real local parallel pool (2 workers) to actually exercise parfor...\n');
pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool('local', 2);

S = load(fullfile(softlubeDir, 'solid_endothelium_P300.mat'));
mesh = S.meshE;
par = S.par;
par.useViscoelasticEndothelium = true;
par.etaE = 1;
par.dt = 1e-5;
par.useObjectiveKelvinVoigt = true;
mesh = prepare_axisym_mesh_cache(mesh);

ndof = size(mesh.nodes,1) * 2;
rng(0);
u    = 1e-15 * randn(ndof, 1);
uOld = 1e-15 * randn(ndof, 1);

tol = 1e-10;
allPass = true;

fprintf('\n--- assemble_finite_def_axisym ---\n');
addpath(backupDir); rehash path;
[Fint_orig, K_orig] = assemble_finite_def_axisym(mesh, u, par);
rmpath(backupDir); rehash path;
[Fint_new, K_new] = assemble_finite_def_axisym(mesh, u, par);
fDiff = max(abs(Fint_orig - Fint_new));
kDiff = max(abs(nonzeros(K_orig - K_new)));
if isempty(kDiff), kDiff = 0; end
fprintf('max |Fint diff| = %.3e, max |K diff| = %.3e\n', fDiff, kDiff);
if fDiff >= tol || kDiff >= tol
    allPass = false;
    fprintf('FAIL\n');
else
    fprintf('PASS\n');
end

fprintf('\n--- assemble_axisym_kelvin_voigt_viscous ---\n');
addpath(backupDir); rehash path;
[Fvisc_orig, Kvisc_orig] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
rmpath(backupDir); rehash path;
[Fvisc_new, Kvisc_new] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
fDiff = max(abs(Fvisc_orig - Fvisc_new));
kNz = nonzeros(Kvisc_orig - Kvisc_new);
kDiff = 0; if ~isempty(kNz), kDiff = max(abs(kNz)); end
fprintf('max |Fvisc diff| = %.3e, max |Kvisc diff| = %.3e\n', fDiff, kDiff);
if fDiff >= tol || kDiff >= tol
    allPass = false;
    fprintf('FAIL\n');
else
    fprintf('PASS\n');
end

fprintf('\n--- assemble_axisym_kelvin_voigt_viscous_force_only ---\n');
addpath(backupDir); rehash path;
Fvisc2_orig = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par);
rmpath(backupDir); rehash path;
Fvisc2_new = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par);
fDiff = max(abs(Fvisc2_orig - Fvisc2_new));
fprintf('max |Fvisc diff| = %.3e\n', fDiff);
if fDiff >= tol
    allPass = false;
    fprintf('FAIL\n');
else
    fprintf('PASS\n');
end

fprintf('\n--- assemble_finite_def_internal_force_only ---\n');
addpath(backupDir); rehash path;
Fint2_orig = assemble_finite_def_internal_force_only(mesh, u, par);
rmpath(backupDir); rehash path;
Fint2_new = assemble_finite_def_internal_force_only(mesh, u, par);
fDiff = max(abs(Fint2_orig - Fint2_new));
fprintf('max |Fint diff| = %.3e\n', fDiff);
if fDiff >= tol
    allPass = false;
    fprintf('FAIL\n');
else
    fprintf('PASS\n');
end

fprintf('\n=== overall: %s ===\n', ternary_local(allPass, 'ALL PASS -- safe to resubmit to the cluster', 'AT LEAST ONE FAILED -- do not resubmit yet'));

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
