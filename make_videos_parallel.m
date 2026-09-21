softlubeDir = '';
if isempty(softlubeDir)
    softlubeDir = fileparts(mfilename('fullpath'));
end
if ~isempty(softlubeDir)
    addpath(softlubeDir);
end

% Ensure all parallel workers have the current directory in their path
pctRunOnAll addpath(softlubeDir);
% Suppress hardware graphics warning across all parallel workers
pctRunOnAll warning('off', 'MATLAB:graphics:noGraphicsAcceleration');

inputname = fullfile(softlubeDir, 'case_7q_restarted.mat');

% Call the parallel generator (file, label, output_dir, fps)
generate_heatmap_videos_parallel(inputname, 'case_7_q_restarted_parallel', 'videos/', 5);
