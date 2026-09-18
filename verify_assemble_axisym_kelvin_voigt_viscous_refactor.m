%% VERIFY_ASSEMBLE_AXISYM_KELVIN_VOIGT_VISCOUS_REFACTOR
% Standalone numerical check: does the parfor-safe rewrite of
% assemble_axisym_kelvin_voigt_viscous.m (accumarray-based Fvisc,
% fixed-offset K slices) produce the same Fvisc/Kvisc as the original, on a
% real mesh/state? Cheap: one assembly call, not a full solve.
%
% Only exercises the cached branch (useCache=true), since that is the only
% branch ever reached in practice -- see the note at the top of
% assemble_axisym_kelvin_voigt_viscous_parfor_safe.m.

clc;
clear all;
softlubeDir = fileparts(mfilename('fullpath'));
addpath(softlubeDir);

fprintf('Building a real endothelium mesh/state to test the viscous assembly refactor on...\n');

S = load(fullfile(softlubeDir, 'solid_endothelium_P300.mat'));
mesh = S.meshE;
par = S.par;

if ~isfield(mesh, 'nelem') || ~isfield(mesh, 'nodes')
    error('meshE does not look like a mesh (no nelem/nodes fields).');
end

% This .mat file's par struct doesn't carry the viscous-specific fields
% (they're set by the run script / cfg at solve time, not stored here) --
% add them so the viscous branch is actually exercised, not short-circuited.
par.useViscoelasticEndothelium = true;
par.etaE = 1;
par.dt = 1e-5;
par.useObjectiveKelvinVoigt = true;

% Build the mesh cache, same as softlube_run_case_global_coupled.m always
% does before assembly -- this is the only branch reached in practice.
mesh = prepare_axisym_mesh_cache(mesh);

ndof = size(mesh.nodes,1) * 2;
rng(0);
% Small enough to avoid any element-inversion issues, nonzero so the
% accumulation/assembly logic is actually exercised.
u    = 1e-15 * randn(ndof, 1);
uOld = 1e-15 * randn(ndof, 1);

fprintf('Running ORIGINAL assemble_axisym_kelvin_voigt_viscous...\n');
[Fvisc_orig, Kvisc_orig] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);

fprintf('Running PARFOR-SAFE assemble_axisym_kelvin_voigt_viscous_parfor_safe...\n');
[Fvisc_new, Kvisc_new] = assemble_axisym_kelvin_voigt_viscous_parfor_safe(mesh, u, uOld, par);

fprintf('\n=== Comparison ===\n');
fDiff = max(abs(Fvisc_orig - Fvisc_new));
fprintf('max |Fvisc_orig - Fvisc_new| = %.6e (scale: max|Fvisc_orig| = %.6e)\n', ...
    fDiff, max(abs(Fvisc_orig)));

kNz = nonzeros(Kvisc_orig - Kvisc_new);
kDiff = 0;
if ~isempty(kNz)
    kDiff = max(abs(kNz));
end
fprintf('max |Kvisc_orig - Kvisc_new| (nonzeros) = %.6e (scale: max|Kvisc_orig| = %.6e)\n', ...
    kDiff, full(max(abs(Kvisc_orig(:)))));

tol = 1e-10;
if fDiff < tol && kDiff < tol
    fprintf('\nPASS: refactor is numerically identical to the original (within %.0e).\n', tol);
else
    fprintf('\nFAIL: refactor differs from the original by more than %.0e -- do not use it yet.\n', tol);
end
