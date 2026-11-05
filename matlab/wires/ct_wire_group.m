function [g, col] = ct_wire_group(w)
%CT_WIRE_GROUP  Display group of a wire and its colour.
p = ct_palette();
if ~w.simulatable, g = 'not simulated'; col = p.light; return; end
switch w.layer_keys{1}
    case 'PVC', g = 'PVC'; col = p.orange;
    case 'XLPE', g = 'XLPE (+PVC jacket)'; col = p.yellow;
    case 'silicone_rubber', g = 'Silicone rubber'; col = p.blue;
    case 'PTFE', g = 'PTFE'; col = p.violet;
    case 'ETFE', g = 'ETFE'; col = p.aqua;
    case 'fluoropolymer_unspecified', g = 'Fluoroplastic (ThermoThin)'; col = p.magenta;
    case 'PVDF', g = 'PVDF'; col = p.green;
    case 'mPPE', g = 'mPPE (EcoWire)'; col = p.red;
    case 'EPR', g = 'EPR'; col = p.grey;
    otherwise, g = 'other'; col = p.light;
end
end
