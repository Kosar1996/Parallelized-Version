% DO NOT DELETE: referenced as the fsolve residual function handle in
% solve_monolithic_two_solids_fsolve_timestep.m, used only when
% par.useTwoSolidSemiAnalyticalJacobian is false (currently unreachable,
% since softlube_prepare_case.m defaults it true and no run script overrides
% it) -- but still referenced by name, so leaving it deleted is a latent bug
% waiting to happen. Verify call sites (grep) before ever removing this again.

function R = monolithic_two_solids_residual_scaled( ...
    y, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
    JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget)
% Residual for the fully coupled two-solid problem.
% The residual is scaled so fsolve sees comparable magnitudes.

    N = numel(z);
    nE = numel(freeE);
    nL = numel(freeL);

    try
        [uE, uL, p] = unpack_two_solid_y( ...
            y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshE, uE, 'endothelium fsolve iterate'), par);
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshL, uL, 'leukocyte fsolve iterate'), par);

        [RE, RL, RF] = monolithic_two_solids_residual_unscaled( ...
            uE, uL, p, old, meshE, interfaceE, baseE, ...
            meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

        R = [
            RE / solidTargetE
            RL / solidTargetL
            RF / fluidTarget
        ];

        if any(~isfinite(R))
            R = 1e12 * ones(nE+nL+N-2,1);
        end
    catch
        % Penalize invalid iterates, especially gap-violating iterates.
        R = 1e12 * ones(nE+nL+N-2,1);
    end
end