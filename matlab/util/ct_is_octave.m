function tf = ct_is_octave()
%CT_IS_OCTAVE  True when running in GNU Octave.
tf = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
