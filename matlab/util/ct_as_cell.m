function c = ct_as_cell(x)
%CT_AS_CELL  jsondecode returns a struct array or a cell array depending on the fields: make it a cell.
if isstruct(x), c = num2cell(x); elseif iscell(x), c = x; else, c = {x}; end
c = c(:)';
end
