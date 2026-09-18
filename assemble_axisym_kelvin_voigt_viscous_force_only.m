function Fvisc = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par)
% Parallelized version (parfor over elements). Profiled as 8.0% of total
% wall-clock time in the 1-timestep single-processor baseline -- one of the
% four hot assembly functions targeted for parallelization. The pre-parfor
% serial version is kept in pre_parallelization_backup/ for reference, and
% the byte-identical verification is in
% verify_assemble_axisym_kelvin_voigt_viscous_force_only_refactor.m.
%
% Fvisc used to be built by direct indexed accumulation (Fvisc(dofs) =
% Fvisc(dofs) + fe), unsafe under parfor because adjacent elements share
% nodes. Fixed by collecting each element's [dofs, fe] into element-exclusive
% slices, summed once via accumarray after the loop.
%
% Note (pre-existing, not introduced here): the non-cached branch below calls
% kelvin_voigt_element_residual_only, which is defined only as a subfunction
% inside softlube_run_case_global_coupled.m -- not visible from this file
% (MATLAB subfunctions are file-local). In practice this is never hit, same
% as the other assembly functions in this folder: softlube_run_case_global_
% coupled.m always builds mesh.axisymCache before assembly runs. Left as
% broken as the original; fixing it is a separate, pre-existing issue.

    ndof = size(mesh.nodes,1)*2;
    Fvisc = zeros(ndof,1);

    if ~(isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium)
        return;
    end
    if ~isfield(par, 'etaE') || par.etaE <= 0
        return;
    end
    useCache = isfield(mesh, 'axisymCache');
    if useCache
        cache = mesh.axisymCache;
    end

    % Cell-array sliced output: see the note in assemble_finite_def_axisym.m
    % for why a computed-range slice (loc = ...; A(loc) = ...) isn't
    % parfor-classifiable and this cell-array indirection is needed instead.
    dofsCell = cell(mesh.nelem, 1);
    feCell = cell(mesh.nelem, 1);

    parfor e = 1:mesh.nelem
        if useCache
            dofs = cache.dofs(e,:).';
            fe = kelvin_voigt_element_residual_only_cached( ...
                cache, e, u(dofs), uOld(dofs), par);
        else
            conn = mesh.conn(e,:);
            Xe   = mesh.nodes(conn,:);
            dofs = reshape([2*conn-1; 2*conn], [], 1);
            fe = kelvin_voigt_element_residual_only(Xe, u(dofs), uOld(dofs), mesh, par);
        end
        dofsCell{e} = dofs;
        feCell{e} = fe;
    end

    iF = zeros(mesh.nelem * 8, 1);
    vF = zeros(mesh.nelem * 8, 1);
    for e = 1:mesh.nelem
        loc = (8*(e-1)+1):(8*e);
        iF(loc) = dofsCell{e};
        vF(loc) = feCell{e};
    end

    Fvisc = accumarray(iF, vF, [ndof, 1]);
end

% Copied unchanged from assemble_axisym_kelvin_voigt_viscous_force_only.m
% (MATLAB subfunctions are only visible within their own file).

function fe = kelvin_voigt_element_residual_only_cached(cache, e, ue, ueOld, par)
    if isfield(par, 'useObjectiveKelvinVoigt') && par.useObjectiveKelvinVoigt
        fe = kelvin_voigt_objective_element_residual_only_cached(cache, e, ue, ueOld, par);
        return;
    end

    fe = zeros(8,1);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnodOld = Rnod + ueOld(1:2:end);
    znodOld = Znod + ueOld(2:2:end);

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

        Fold = deformation_gradient_from_nodal([], rnodOld, znodOld, N, dNdX, Rg);
        F = current_deformation_gradient_from_cached(cache, e, ue, g);
        Pvisc = par.etaE * (F - Fold) / par.dt;
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

            fe(2*a) = fe(2*a) + ...
                (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
        end
    end
end

function fe = kelvin_voigt_objective_element_residual_only_cached(cache, e, ue, ueOld, par)
    fe = zeros(8,1);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnodOld = Rnod + ueOld(1:2:end);
    znodOld = Znod + ueOld(2:2:end);

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

        Fold = deformation_gradient_from_nodal([], rnodOld, znodOld, N, dNdX, Rg);
        F = current_deformation_gradient_from_cached(cache, e, ue, g);

        Pvisc = objective_kelvin_voigt_piola(F, Fold, par);
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

            fe(2*a) = fe(2*a) + ...
                (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
        end
    end
end
