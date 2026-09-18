%% VERIFY_ASSEMBLE_AXISYM_KELVIN_VOIGT_VISCOUS_FORCE_ONLY_REFACTOR
% Standalone numerical check: does the parfor-safe rewrite of
% assemble_axisym_kelvin_voigt_viscous_force_only.m (accumarray-based
% Fvisc) produce the same Fvisc as the original, on a real mesh/state?
% Cheap: one assembly call, not a full solve.

clc;
clear all;
softlubeDir = fileparts(mfilename('fullpath'));
addpath(softlubeDir);

fprintf('Building a real endothelium mesh/state to test the viscous force-only refactor on...\n');

S = load(fullfile(softlubeDir, 'solid_endothelium_P300.mat'));
mesh = S.meshE;
par = S.par;

if ~isfield(mesh, 'nelem') || ~isfield(mesh, 'nodes')
    error('meshE does not look like a mesh (no nelem/nodes fields).');
end

par.useViscoelasticEndothelium = true;
par.etaE = 1;
par.dt = 1e-5;
par.useObjectiveKelvinVoigt = true;

mesh = prepare_axisym_mesh_cache(mesh);

ndof = size(mesh.nodes,1) * 2;
rng(0);
u    = 1e-15 * randn(ndof, 1);
uOld = 1e-15 * randn(ndof, 1);

fprintf('Running ORIGINAL assemble_axisym_kelvin_voigt_viscous_force_only...\n');
Fvisc_orig = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par);

fprintf('Running PARFOR-SAFE assemble_axisym_kelvin_voigt_viscous_force_only_parfor_safe...\n');
Fvisc_new = assemble_axisym_kelvin_voigt_viscous_force_only_parfor_safe(mesh, u, uOld, par);

fprintf('\n=== Comparison ===\n');
fDiff = max(abs(Fvisc_orig - Fvisc_new));
fprintf('max |Fvisc_orig - Fvisc_new| = %.6e (scale: max|Fvisc_orig| = %.6e)\n', ...
    fDiff, max(abs(Fvisc_orig)));

tol = 1e-10;
if fDiff < tol
    fprintf('\nPASS: refactor is numerically identical to the original (within %.0e).\n', tol);
else
    fprintf('\nFAIL: refactor differs from the original by more than %.0e -- do not use it yet.\n', tol);
end
