%% VERIFY_ASSEMBLE_FINITE_DEF_INTERNAL_FORCE_ONLY_REFACTOR
% Standalone numerical check: does the parfor-safe rewrite of
% assemble_finite_def_internal_force_only.m (accumarray-based Fint) produce
% the same Fint as the original, on a real mesh/state? Cheap: one assembly
% call, not a full solve.

clc;
clear all;
softlubeDir = fileparts(mfilename('fullpath'));
addpath(softlubeDir);

fprintf('Building a real endothelium mesh/state to test the internal-force refactor on...\n');

S = load(fullfile(softlubeDir, 'solid_endothelium_P300.mat'));
mesh = S.meshE;
par = S.par;

if ~isfield(mesh, 'nelem') || ~isfield(mesh, 'nodes')
    error('meshE does not look like a mesh (no nelem/nodes fields).');
end
if ~isfield(par, 'Ge') || ~isfield(par, 'Ke')
    error('par (from the .mat file) is missing Ge/Ke -- needed by the element physics.');
end

% Build the mesh cache, same as softlube_run_case_global_coupled.m always
% does before assembly -- this is the only branch reached in practice.
mesh = prepare_axisym_mesh_cache(mesh);

ndof = size(mesh.nodes,1) * 2;
rng(0);
u = 1e-15 * randn(ndof, 1);

fprintf('Running ORIGINAL assemble_finite_def_internal_force_only...\n');
Fint_orig = assemble_finite_def_internal_force_only(mesh, u, par);

fprintf('Running PARFOR-SAFE assemble_finite_def_internal_force_only_parfor_safe...\n');
Fint_new = assemble_finite_def_internal_force_only_parfor_safe(mesh, u, par);

fprintf('\n=== Comparison ===\n');
fDiff = max(abs(Fint_orig - Fint_new));
fprintf('max |Fint_orig - Fint_new| = %.6e (scale: max|Fint_orig| = %.6e)\n', ...
    fDiff, max(abs(Fint_orig)));

tol = 1e-10;
if fDiff < tol
    fprintf('\nPASS: refactor is numerically identical to the original (within %.0e).\n', tol);
else
    fprintf('\nFAIL: refactor differs from the original by more than %.0e -- do not use it yet.\n', tol);
end
