function root = setup_paths()
%SETUP_PATHS  Add the cable-temperature MATLAB/Octave toolbox folders to the path.
%   Run once per session (all scripts call it).  Works in MATLAB R2016b+ and GNU Octave 7+.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'core'), fullfile(root, 'heat'), fullfile(root, 'iec'), ...
        fullfile(root, 'wires'), fullfile(root, 'util'), fullfile(root, 'simulink'));
end
