%% VERIFY_ASSEMBLE_FINITE_DEF_AXISYM_REFACTOR
% Standalone numerical check: does the parfor-safe rewrite of
% assemble_finite_def_axisym.m (accumarray-based Fint, fixed-offset K
% slices) produce the same Fint/K as the original, on a real mesh/state?
%
% Run this BEFORE trusting assemble_finite_def_axisym_parfor_safe.m enough to
% replace the live file. Cheap: one assembly call, not a full solve. Both
% functions already sit side by side in this folder under their own names,
% so no setup is needed -- just run this script.

clc;
clear all;
softlubeDir = fileparts(mfilename('fullpath'));
addpath(softlubeDir);

fprintf('Building a real endothelium mesh/state to test the assembly refactor on...\n');

S = load(fullfile(softlubeDir, 'solid_endothelium_P300.mat'));
mesh = S.meshE;
par = S.par;

if ~isfield(mesh, 'nelem') || ~isfield(mesh, 'nodes')
    error(['meshE does not look like a mesh (no nelem/nodes fields). ', ...
        'Inspect solid_endothelium_P300.mat manually and point this script at the ', ...
        'right field before relying on this check.']);
end
if ~isfield(par, 'Ge') || ~isfield(par, 'Ke')
    error('par (from the .mat file) is missing Ge/Ke -- needed by the element physics.');
end

% Build the mesh cache, same as softlube_run_case_global_coupled.m always
% does before assembly -- this is the only branch reached in practice.
mesh = prepare_axisym_mesh_cache(mesh);

ndof = size(mesh.nodes,1) * 2;
rng(0);
% Small enough to guarantee no element inversion regardless of local mesh
% scale, while still nonzero (this test only needs to exercise the assembly
% and accumulation logic, not represent anything physically meaningful).
u = 1e-15 * randn(ndof, 1);

fprintf('Running ORIGINAL assemble_finite_def_axisym...\n');
[Fint_orig, K_orig] = assemble_finite_def_axisym(mesh, u, par);

fprintf('Running PARFOR-SAFE assemble_finite_def_axisym_parfor_safe...\n');
[Fint_new, K_new] = assemble_finite_def_axisym_parfor_safe(mesh, u, par);

fprintf('\n=== Comparison ===\n');
fDiff = max(abs(Fint_orig - Fint_new));
fprintf('max |Fint_orig - Fint_new| = %.6e (Fint scale: max|Fint_orig| = %.6e)\n', ...
    fDiff, max(abs(Fint_orig)));

KDiff = max(abs(nonzeros(K_orig - K_new)));
if isempty(KDiff)
    KDiff = 0;
end
fprintf('max |K_orig - K_new| (nonzeros)   = %.6e (K scale: max|K_orig| = %.6e)\n', ...
    KDiff, full(max(abs(K_orig(:)))));

tol = 1e-10;
if fDiff < tol && KDiff < tol
    fprintf('\nPASS: refactor is numerically identical to the original (within %.0e).\n', tol);
else
    fprintf('\nFAIL: refactor differs from the original by more than %.0e -- do not use it yet.\n', tol);
end
