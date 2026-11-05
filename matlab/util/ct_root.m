function root = ct_root()
%CT_ROOT  Folder of the MATLAB toolbox (parent of util/).
root = fileparts(fileparts(mfilename('fullpath')));
end
