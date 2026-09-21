% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Reference Frame Divergence Fix:
%    Preserved uOld across all load increments so Kelvin-Voigt viscous forces 
%    evaluate against the true global time-step baseline.
%
% 2. Integrated Viscous Assembly:
%    Directly passes uOld to assemble_finite_def_axisym.m to evaluate total
%    Piola stress (P_elastic + P_visc) without double-assembly calls.
%
% 3. Line-Search Scale Integration:
%    Passes parTrial with parTrial.alpha_ls into internal force assembly
%    during line-search trial steps.
%
% 4. Flexible 8-Argument Signature Fix (ISOLATED TEST):
%    Updated function header to accept uInitialOrMesh, preventing 
%    "Too many input arguments" errors when callers pass an 8th argument.
% -------------------------------------------------------------------------

function uNew = solve_finite_def_solid(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitialOrMesh)
    nLoadSteps = 5;
    if isfield(par, 'solidLoadSteps') && isfinite(par.solidLoadSteps) && par.solidLoadSteps >= 1
        nLoadSteps = max(1, round(par.solidLoadSteps));
    end

    uInitial = [];
    if nargin >= 8 && ~isempty(uInitialOrMesh)
        if isnumeric(uInitialOrMesh)
            uInitial = uInitialOrMesh;
        end
    end

    if isempty(uInitial)
        uCurrent = uOld;
    else
        uCurrent = uInitial;
    end

    if nLoadSteps <= 1
        uNew = solve_step_with_fallback(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uCurrent);
        return;
    end

    for step = 1:nLoadSteps
        frac = step / nLoadSteps;
        tractionStep = traction;
        if isfield(tractionStep, 'normal')
            tractionStep.normal = frac * tractionStep.normal;
        end
        if isfield(tractionStep, 'tangent')
            tractionStep.tangent = frac * tractionStep.tangent;
        end
        uCurrent = solve_step_with_fallback(mesh, uOld, tractionStep, interfaceNodes, baseNodes, supportType, par, uCurrent);
    end
    uNew = uCurrent;
end

function uNew = solve_step_with_fallback(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
    try
        uNew = solve_finite_def_solid_singlestep(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial);
    catch ME
        if ~(contains(ME.message, 'stalled at the trust-region floor') || ...
             contains(ME.message, 'line search failed before equilibrium') || ...
             contains(ME.message, 'hit max iterations before equilibrium') || ...
             contains(ME.message, 'top-of-iteration evaluation hit an inverted element'))
            rethrow(ME);
        end
        uNew = solve_finite_def_solid_fsolve(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial);
    end
end

function uNew = solve_finite_def_solid_fsolve(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
    ndof = size(mesh.nodes,1)*2;
    if nargin >= 8 && ~isempty(uInitial) && isnumeric(uInitial)
        u0 = uInitial;
    else
        u0 = uOld;
    end
    [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u0(fixDofs) = fixVals;

    y0 = u0(free);
    fun = @(y) finite_def_solid_residual_jac(y, uOld, fixDofs, fixVals, free, ndof, ...
        mesh, interfaceNodes, traction, par);

    fsolveMaxIterSolid = 200;
    if isfield(par, 'solidFsolveMaxIterations') && isfinite(par.solidFsolveMaxIterations)
        fsolveMaxIterSolid = par.solidFsolveMaxIterations;
    end
    fsolveMaxFunEvalSolid = 2000;
    if isfield(par, 'solidFsolveMaxFunctionEvaluations') && isfinite(par.solidFsolveMaxFunctionEvaluations)
        fsolveMaxFunEvalSolid = par.solidFsolveMaxFunctionEvaluations;
    end

    opts = optimoptions('fsolve', ...
        'Algorithm', 'trust-region-dogleg', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-8, ...
        'MaxIterations', fsolveMaxIterSolid, ...
        'MaxFunctionEvaluations', fsolveMaxFunEvalSolid);

    [ySol, ~, exitflag, output] = fsolve(fun, y0, opts);

    if isfield(par, 'debugVerbose') && par.debugVerbose
        fprintf('      solid fsolve fallback: exitflag=%d, iterations=%d\n', exitflag, output.iterations);
    end

    if exitflag <= 0
        error(['Finite-deformation solid solve failed in both the custom Newton loop ', ...
               'and the fsolve fallback (fsolve exitflag=%d).'], exitflag);
    end

    uNew = uOld;
    uNew(free) = ySol;
    uNew(fixDofs) = fixVals;
end

function [R, J] = finite_def_solid_residual_jac(y, uOld, fixDofs, fixVals, free, ndof, ...
    mesh, interfaceNodes, traction, par)
u = uOld;
u(free) = y;
u(fixDofs) = fixVals;

Fext = zeros(ndof,1);
try
    [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
    [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par, uOld);
catch ME
    if ~(contains(ME.message, 'Negative or zero J') || ...
         contains(ME.message, 'Non-positive radius') || ...
         contains(ME.message, 'Element inverted'))
        rethrow(ME);
    end
    n = numel(free);
    R = 1e6 * ones(n,1);
    if nargout > 1
        J = speye(n) * 1e6;
    end
    return;
end

refNormFloor = 1e-6;
if isfield(par, 'solidRefNormFloor') && isfinite(par.solidRefNormFloor)
    refNormFloor = par.solidRefNormFloor;
end
refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), refNormFloor]);
R = (Fint(free) - Fext(free)) / refNorm;
if nargout > 1
    Ktot = Ktan - Kext;
    J = Ktot(free, free) / refNorm;
end
end

function uNew = solve_finite_def_solid_singlestep(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
    ndof = size(mesh.nodes,1)*2;
    if nargin >= 8 && ~isempty(uInitial) && isnumeric(uInitial)
        u = uInitial;
    else
        u = uOld;
    end

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    absTol = 1e-13;
    if isfield(par, 'solidAbsTol')
        absTol = par.solidAbsTol;
    end

    refNormFloor = 1e-6;
    if isfield(par, 'solidRefNormFloor') && isfinite(par.solidRefNormFloor)
        refNormFloor = par.solidRefNormFloor;
    end

    fallbackAbsTol = [];
    if isfield(par, 'solidFallbackAbsTol')
        fallbackAbsTol = par.solidFallbackAbsTol;
    end
    useFallbackAbsTol = ~isempty(fallbackAbsTol) && ...
        isfinite(fallbackAbsTol) && fallbackAbsTol > 0;

    fallbackRelTol = inf;
    if isfield(par, 'solidFallbackRelTol')
        fallbackRelTol = par.solidFallbackRelTol;
    end
    fallbackMinIterations = 2;
    if isfield(par, 'solidFallbackMinIterations')
        fallbackMinIterations = max(1, round(par.solidFallbackMinIterations));
    end

    trustU = par.trustU0;
    if isfield(par, 'solidTrustU0')
        trustU = par.solidTrustU0;
    end

    trustUMin = par.trustUMin;
    if isfield(par, 'solidTrustUMin')
        trustUMin = par.solidTrustUMin;
    end

    trustUMax = par.trustUMax;
    if isfield(par, 'solidTrustUMax')
        trustUMax = par.solidTrustUMax;
    end

    bestRel = inf;
    bestAbs = inf;
    bestU = u;
    bestUpdated = false;

    stallFloorCount = 0;
    stallRelAtFloorStart = inf;
    stallMaxIters = 15;

    for it = 1:par.newtonMaxItSolid

        try
            Fext = zeros(ndof,1);
            [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
            [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par, uOld);
        catch ME
            if contains(ME.message, 'Negative or zero J') || ...
               contains(ME.message, 'Non-positive radius') || ...
               contains(ME.message, 'Element inverted')
                error(['Solid Newton top-of-iteration evaluation hit an inverted element ', ...
                       '(iteration %d): %s'], it, ME.message);
            else
                rethrow(ME);
            end
        end

        R = Fint - Fext;
        Ktot = Ktan - Kext;

        Rf  = R(free);
        Kff = Ktot(free, free);

        resNorm = norm(Rf, inf);
        refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), refNormFloor]);
        relNorm = resNorm / refNorm;

        if isfield(par, 'debugVerbose') && par.debugVerbose
            fprintf('      solid Newton %d: rel=%.3e abs=%.3e trustU=%.3e\n', ...
                it, relNorm, resNorm, trustU);
        end

        if relNorm < bestRel
            bestRel = relNorm;
            bestAbs = resNorm;
            bestU = u;
        end

        if relNorm < par.newtonTolSolid || resNorm < absTol
            uNew = u;
            return;
        end

        if useFallbackAbsTol && it >= fallbackMinIterations && ...
                resNorm < fallbackAbsTol && relNorm < fallbackRelTol
            uNew = u;
            return;
        end

        dscale = sqrt(abs(full(diag(Kff))));
        dscale(dscale < eps(class(dscale)) | ~isfinite(dscale)) = 1;
        Dinv = spdiags(1./dscale, 0, numel(dscale), numel(dscale));
        Kscaled = Dinv * Kff * Dinv;
        Rscaled = Dinv * Rf;
        du_scaled = -Kscaled \ Rscaled;
        du_free = dscale .\ du_scaled;
        if ~all(isfinite(du_free))
            error('Solid Newton produced non-finite displacement increments.');
        end

        duMax = max(abs(du_free));
        if duMax > trustU
            stepScale = trustU / duMax;
        else
            stepScale = 1.0;
        end

        alpha = 1.0;
        accepted = false;
        resAccepted = inf;

        for ls = 1:par.lineSearchMax
            uTrial = u;
            uTrial(free) = uTrial(free) + alpha * stepScale * du_free;
            uTrial(fixDofs) = fixVals;

            parTrial = par;
            parTrial.alpha_ls = alpha * stepScale;

            try
                FextTrial = zeros(ndof,1);
                [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);
                [FintTrial, ~] = assemble_finite_def_axisym(mesh, uTrial, parTrial, uOld);
                
                Rtrial = FintTrial - FextTrial;
                resTrial = norm(Rtrial(free), inf);
                refTrial = max([norm(FextTrial(free), inf), norm(FintTrial(free), inf), refNormFloor]);
                relTrial = resTrial / refTrial;

                if relTrial < bestRel
                    bestRel = relTrial;
                    bestAbs = resTrial;
                    bestU = uTrial;
                    bestUpdated = true;
                end

                if resTrial < resNorm || resTrial < absTol
                    u = uTrial;
                    resAccepted = resTrial;
                    accepted = true;
                    break;
                end

            catch ME
                if contains(ME.message, 'Negative or zero J') || ...
                   contains(ME.message, 'Non-positive radius') || ...
                   contains(ME.message, 'Element inverted')
                else
                    rethrow(ME);
                end
            end

            alpha = 0.5 * alpha;
        end

        if ~accepted
            trustU = 0.5 * trustU;
            if trustU < trustUMin
                error(['Solid Newton line search failed before equilibrium. ', ...
                       'best relative residual = %.3e, best absolute residual = %.3e.'], ...
                       bestRel, bestAbs);
            end
            continue;
        end

        if resAccepted < 0.25 * resNorm
            trustU = min(2.0 * trustU, trustUMax);
        elseif alpha < 0.25 || stepScale < 0.25
            trustU = max(0.5 * trustU, trustUMin);
        end

        if trustU <= 1.001 * trustUMin
            if stallFloorCount == 0
                stallRelAtFloorStart = bestRel;
            end
            stallFloorCount = stallFloorCount + 1;
            relImprovement = (stallRelAtFloorStart - bestRel) / max(stallRelAtFloorStart, 1e-30);
            if stallFloorCount >= stallMaxIters && relImprovement < 1e-3
                error(['Finite-deformation solid solve stalled at the trust-region floor: ', ...
                       '%d iterations with no meaningful progress (best relative residual = %.3e, ', ...
                       'best absolute residual = %.3e). Bailed out early instead of exhausting ', ...
                       'newtonMaxItSolid.'], stallFloorCount, bestRel, bestAbs);
            end
        else
            stallFloorCount = 0;
            stallRelAtFloorStart = inf;
        end
    end

    if useFallbackAbsTol && bestUpdated && ...
            bestAbs < fallbackAbsTol && bestRel < fallbackRelTol
        uNew = bestU;
        return;
    end

    error(['Finite-deformation solid solve hit max iterations before equilibrium. ', ...
           'best relative residual = %.3e, best absolute residual = %.3e.'], ...
           bestRel, bestAbs);
end