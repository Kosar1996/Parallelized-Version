function [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par)
% Parallelized version (parfor over elements). Profiled as 54.2% of total
% wall-clock time in the 1-timestep single-processor baseline -- the single
% hottest function, and the top target for parallelization. The pre-parfor
% serial version is kept in pre_parallelization_backup/ for reference, and
% the byte-identical verification is in
% verify_assemble_axisym_kelvin_voigt_viscous_refactor.m.
%
% Same fix as assemble_finite_def_axisym.m: Fvisc used to be built by direct
% indexed accumulation (Fvisc(dofs) = Fvisc(dofs) + fe), unsafe under parfor
% because adjacent elements share nodes; and the non-cached branch used a
% running ptr counter, a loop-carried dependency, also unsafe under parfor.
% Both fixed the same way -- element-exclusive slices, summed once via
% accumarray after the loop.
%
% Note (pre-existing, not introduced here): the non-cached branch below calls
% kelvin_voigt_element_residual_tangent, which is not defined anywhere in
% this codebase. In practice this is never hit -- softlube_run_case_global_
% coupled.m always builds mesh.axisymCache before assembly runs -- but this
% branch is left exactly as broken as the original, since fixing it is a
% separate, pre-existing issue and out of scope for parallelization prep.

    ndof = size(mesh.nodes,1)*2;
    Fvisc = zeros(ndof,1);

    if ~(isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium)
        Kvisc = sparse(ndof, ndof);
        return;
    end
    if ~isfield(par, 'etaE') || par.etaE <= 0
        Kvisc = sparse(ndof, ndof);
        return;
    end

    useCache = isfield(mesh, 'axisymCache');
    if useCache
        cache = mesh.axisymCache;
        iK = cache.iK;
        jK = cache.jK;
    else
        nnzLocal = mesh.nelem * 64;
        iK = zeros(nnzLocal,1);
        jK = zeros(nnzLocal,1);
    end
    vK = zeros(size(iK));

    % Cell-array sliced output: see the note in assemble_finite_def_axisym.m
    % for why a computed-range slice (loc = ...; A(loc) = ...) isn't
    % parfor-classifiable and this cell-array indirection is needed instead.
    dofsCell = cell(mesh.nelem, 1);
    feCell = cell(mesh.nelem, 1);
    KeCell = cell(mesh.nelem, 1);
    if ~useCache
        iiCell = cell(mesh.nelem, 1);
        jjCell = cell(mesh.nelem, 1);
    end

    parfor e = 1:mesh.nelem
        if useCache
            dofs = cache.dofs(e,:).';
            [fe, Ke] = kelvin_voigt_element_residual_tangent_cached( ...
                cache, e, u(dofs), uOld(dofs), par);
        else
            conn = mesh.conn(e,:);
            Xe   = mesh.nodes(conn,:);
            dofs = reshape([2*conn-1; 2*conn], [], 1);
            [fe, Ke] = kelvin_voigt_element_residual_tangent( ...
                Xe, u(dofs), uOld(dofs), mesh, par);
        end

        dofsCell{e} = dofs;
        feCell{e} = fe;
        KeCell{e} = Ke;
        if ~useCache
            [ii, jj] = ndgrid(dofs, dofs);
            iiCell{e} = ii(:);
            jjCell{e} = jj(:);
        end
    end

    iF = zeros(mesh.nelem * 8, 1);
    vF = zeros(mesh.nelem * 8, 1);

    for e = 1:mesh.nelem
        locF = (8*(e-1)+1):(8*e);
        iF(locF) = dofsCell{e};
        vF(locF) = feCell{e};

        locK = (64*(e-1)+1):(64*e);
        if ~useCache
            iK(locK) = iiCell{e};
            jK(locK) = jjCell{e};
        end
        vK(locK) = KeCell{e}(:);
    end

    Fvisc = accumarray(iF, vF, [ndof, 1]);
    Kvisc = sparse(iK, jK, vK, ndof, ndof);
end

% Copied unchanged from assemble_axisym_kelvin_voigt_viscous.m (MATLAB
% subfunctions are only visible within their own file).

function [fe, Ke] = kelvin_voigt_element_residual_tangent_cached(cache, e, ue, ueOld, par)
    fe = zeros(8,1);
    Ke = zeros(8,8);

    if isfield(par, 'useObjectiveKelvinVoigt') && par.useObjectiveKelvinVoigt
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
            [Pvisc, kv] = objective_kelvin_voigt_piola(F, Fold, par);
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

            for alpha = 1:8
                dP = objective_kelvin_voigt_piola_tangent( ...
                    local_dF_from_dof(alpha, N, dNdX, Rg), kv, par);

                for a = 1:4
                    dNa_dR = dNdX(a,1);
                    dNa_dZ = dNdX(a,2);
                    Na     = N(a);

                    Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                        (dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg)) * Wgp;

                    Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                        (dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ) * Wgp;
                end
            end
        end
        return;
    end

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

        for alpha = 1:8
            dP = (par.etaE/par.dt) * local_dF_from_dof(alpha, N, dNdX, Rg);

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    (dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg)) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    (dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ) * Wgp;
            end
        end
    end
end
