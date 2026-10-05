%% ANALYZE_0930 -- builds the numbers for Master Tables 2, 3, 4
% Run from 0930_parallel/ after all jobs finish (0930_serial/ must sit next
% to it). Reads:
%   ../0930_serial/out_0930_serial_base.mat        serial baseline (original code)
%   ./out_0930_par_N{1,2,4,8,16}.mat               parallel runs
%   ../0930_serial/profile_data_0930_serial_prof.mat   (Table 4, serial)
%   ./profile_data_0930_par_prof_N8.mat                 (Table 4, parallel)
% Writes (into ./analysis_0930/):
%   table2_correctness.csv, table3_scaling.csv, table4_profiling.csv,
%   step_times.csv, analysis_0930_log.txt

clc; clear all;
here = fileparts(mfilename('fullpath'));
serialDir = fullfile(here, '..', '0930_serial');
outDir = fullfile(here, 'analysis_0930');
if ~exist(outDir, 'dir'), mkdir(outDir); end
diary(fullfile(outDir, 'analysis_0930_log.txt'));

runs = {'serial_base', fullfile(serialDir, 'out_0930_serial_base.mat'), 0};
for n = [1 2 4 8 16]
    runs(end+1,:) = {sprintf('par_N%d', n), fullfile(here, sprintf('out_0930_par_N%d.mat', n)), n}; %#ok<SAGROW>
end

R = struct([]); r = struct('tag',[],'N',[],'out',[],'wall',[],'fromCheckpoint',[]);
for i = 1:size(runs,1)
    f = runs{i,2};
    if ~exist(f, 'file')
        fprintf('MISSING: %s\n', f); continue;
    end
    d = load(f);
    r.tag = runs{i,1}; r.N = runs{i,3};
    r.out = d.out;
    % Normal finish: wallClockSeconds saved by run_0930_bench.m.
    % Time-limit kill: only the per-step checkpoint exists (no wallClockSeconds);
    % then use the sum of per-step wall times and flag it.
    if isfield(d, 'wallClockSeconds')
        r.wall = d.wallClockSeconds; r.fromCheckpoint = false;
    else
        r.wall = sum(d.out.stepWallTimeHist(:), 'omitnan'); r.fromCheckpoint = true;
        fprintf('NOTE: %s has no wallClockSeconds (checkpoint only, stopStep=%d); wall = sum(stepWallTimeHist) = %.1f s\n', ...
            runs{i,1}, d.out.stopStep, r.wall);
    end
    if isempty(R), R = r; else, R(end+1) = r; end %#ok<SAGROW>
end
tags = {R.tag};
iBase = find(strcmp(tags, 'serial_base'));
iN1   = find(strcmp(tags, 'par_N1'));

%% ---------------- Table 2: correctness ----------------
fid = fopen(fullfile(outDir, 'table2_correctness.csv'), 'w');
fprintf(fid, 'run,N,stopStep,final_t,t_array_matches_serial,t_array_matches_N1,maxP2D_final,maxP2D_hist_diff_vs_serial,maxP2D_hist_diff_vs_N1,PHist_full_maxabsdiff_vs_serial,uE_uL_final_maxabsdiff_vs_serial\n');
for i = 1:numel(R)
    o = R(i).out;
    [tS, pS, PHs, uS] = cmp(o, R, iBase);
    [tN, pN] = cmp(o, R, iN1);
    fprintf(fid, '%s,%d,%d,%.9e,%s,%s,%.9e,%s,%s,%s,%s\n', R(i).tag, R(i).N, o.stopStep, o.t(end), ...
        tS, tN, lastfinite(o.p2DMaxHist), pS, pN, PHs, uS);
    fprintf('%-12s stopStep=%d t_end=%.6e | vs serial: t %s, max|P2D| diff %s, PHist %s, uE/uL %s | vs N1: t %s, max|P2D| diff %s\n', ...
        R(i).tag, o.stopStep, o.t(end), tS, pS, PHs, uS, tN, pN);
end
fclose(fid);

%% ---------------- Table 3: scaling ----------------
% Common step range: every run is timed over the SAME steps (1..nC), from its
% per-step wall-clock history, so runs that stopped at different steps
% (time limit) are compared on identical work.
nC = min(arrayfun(@(x) numel(x.out.stepWallTimeHist), R));
for i = 1:numel(R)
    w = R(i).out.stepWallTimeHist(:);
    R(i).wall = sum(w(1:nC), 'omitnan');
end
fprintf('Table 3 timed over common steps 1-%d\n', nC);
fid = fopen(fullfile(outDir, 'table3_scaling.csv'), 'w');
fprintf(fid, 'run,N,wall_s,wall_hr,speedup_vs_serial,speedup_vs_N1,stopStep,avg_s_per_step_excl_last,from_checkpoint\n');
for i = 1:numel(R)
    o = R(i).out;
    spS = R(iBase).wall / R(i).wall;
    spN = NaN; if ~isempty(iN1), spN = R(iN1).wall / R(i).wall; end
    st = NaN;
    if isfield(o, 'stepWallTimeHist') && numel(o.stepWallTimeHist) > 1
        w = o.stepWallTimeHist(:); st = mean(w(1:nC), 'omitnan');
    end
    fprintf(fid, '%s,%d,%.1f,%.3f,%.3f,%.3f,%d,%.1f,%d\n', R(i).tag, R(i).N, R(i).wall, R(i).wall/3600, spS, spN, o.stopStep, st, R(i).fromCheckpoint);
end
fclose(fid);

% per-step wall time (all runs, for the appendix)
fid = fopen(fullfile(outDir, 'step_times.csv'), 'w');
fprintf(fid, 'run,step,t,dt_s_wall\n');
for i = 1:numel(R)
    o = R(i).out;
    if ~isfield(o, 'stepWallTimeHist'), continue; end
    for k = 1:numel(o.stepWallTimeHist)
        fprintf(fid, '%s,%d,%.9e,%.2f\n', R(i).tag, k, o.t(min(k,end)), o.stepWallTimeHist(k));
    end
end
fclose(fid);

%% ---------------- Table 4: profiling ----------------
% Matched 3-step profiling pair (serial vs parallel N=8), used because the
% full 14-step serial profiling job was cancelled to free cluster resources.
profFiles = {fullfile(serialDir, 'profile_data_0930_serial_prof3.mat'), ...
             fullfile(here, 'profile_data_0930_par_prof3_N8.mat')};
T = cell(1,2); W = [NaN NaN];
for j = 1:2
    if ~exist(profFiles{j}, 'file'), fprintf('MISSING: %s\n', profFiles{j}); continue; end
    d = load(profFiles{j}); W(j) = d.wallClockSeconds;
    T{j} = selftimes(d.pInfo);
end
if ~isempty(T{1}) && ~isempty(T{2})
    names = union(T{1}.name, T{2}.name, 'stable');
    rows = {};
    for k = 1:numel(names)
        [c1, s1, t1] = lookup(T{1}, names{k});
        [c2, s2, t2] = lookup(T{2}, names{k});
        rows(end+1,:) = {names{k}, c1, c2, s1, s2, t1, t2}; %#ok<SAGROW>
    end
    selfMax = max(cell2mat(rows(:,4)), cell2mat(rows(:,5)));
    [~, ord] = sort(selfMax, 'descend');
    ord = ord(1:min(20, numel(ord)));
    fid = fopen(fullfile(outDir, 'table4_profiling.csv'), 'w');
    fprintf(fid, 'function,calls_serial,calls_parallel,self_s_serial,self_s_parallel,total_s_serial,total_s_parallel\n');
    fprintf(fid, 'WALL_CLOCK_TOTAL,,,%.1f,%.1f,,\n', W(1), W(2));
    for k = ord(:).'
        fprintf(fid, '"%s",%d,%d,%.1f,%.1f,%.1f,%.1f\n', rows{k,:});
    end
    fclose(fid);
    fprintf('\nTable 4 written (top %d functions by self time).\n', numel(ord));
end
diary off;
fprintf('\nAll CSVs in %s\n', outDir);

%% ================= helpers =================
function [tStr, pStr, PHstr, uStr] = cmp(o, R, iRef)
    % Compares over the COMMON step range (runs may stop at different steps,
    % e.g. N=16 saved 14 steps, the others 13 at the time limit).
    tStr = 'n/a'; pStr = 'n/a'; PHstr = 'n/a'; uStr = 'n/a';
    if isempty(iRef), return; end
    r = R(iRef).out;
    n = min(numel(o.t), numel(r.t));
    if isequal(o.t(1:n), r.t(1:n)), tStr = sprintf('YES exact (steps 1-%d)', n); else, tStr = sprintf('NO (steps 1-%d)', n); end
    if isfield(o,'p2DMaxHist') && isfield(r,'p2DMaxHist')
        m = min(n, min(numel(o.p2DMaxHist), numel(r.p2DMaxHist)));
        pStr = sprintf('%.3e (steps 1-%d)', max(abs(o.p2DMaxHist(1:m) - r.p2DMaxHist(1:m)), [], 'omitnan'), m);
    end
    if isfield(o,'PHist') && isfield(r,'PHist') && size(o.PHist,1)==size(r.PHist,1) && size(o.PHist,2)==size(r.PHist,2)
        m = min(n, min(size(o.PHist,3), size(r.PHist,3)));
        A = o.PHist(:,:,1:m); B = r.PHist(:,:,1:m);
        PHstr = sprintf('%.3e (steps 1-%d)', max(abs(A(:) - B(:)), [], 'omitnan'), m);
    end
    if isfield(o,'stateHist') && isfield(r,'stateHist') && numel(o.stateHist) >= n && numel(r.stateHist) >= n
        du = 0;
        try
            du = max(abs(o.stateHist{n}.uE(:) - r.stateHist{n}.uE(:)));
            if isfield(o.stateHist{n},'uL') && ~isempty(o.stateHist{n}.uL)
                du = max(du, max(abs(o.stateHist{n}.uL(:) - r.stateHist{n}.uL(:))));
            end
            uStr = sprintf('%.3e (step %d)', du, n);
        catch
        end
    end
end

function v = lastfinite(x)
    x = x(isfinite(x)); if isempty(x), v = NaN; else, v = x(end); end
end

function T = selftimes(pInfo)
    FT = pInfo.FunctionTable;
    n = numel(FT);
    T.name = cell(n,1); T.calls = zeros(n,1); T.self = zeros(n,1); T.total = zeros(n,1);
    for k = 1:n
        ch = 0;
        if ~isempty(FT(k).Children), ch = sum([FT(k).Children.TotalTime]); end
        T.name{k}  = FT(k).FunctionName;
        T.calls(k) = FT(k).NumCalls;
        T.total(k) = FT(k).TotalTime;
        T.self(k)  = FT(k).TotalTime - ch;
    end
end

function [c, s, t] = lookup(T, name)
    c = 0; s = 0; t = 0;
    if isempty(T), return; end
    m = strcmp(T.name, name);
    if any(m), c = sum(T.calls(m)); s = sum(T.self(m)); t = sum(T.total(m)); end
end
