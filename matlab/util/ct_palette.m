function p = ct_palette()
%CT_PALETTE  Colours used by the figures (same as the Python version).
h = @(s) [hex2dec(s(2:3)) hex2dec(s(4:5)) hex2dec(s(6:7))] / 255;
p.orange = h('#eb6834'); p.blue = h('#2a78d6'); p.grey = h('#6b6a65'); p.violet = h('#6250d6');
p.aqua = h('#1baf7a'); p.red = h('#e34948'); p.yellow = h('#eda100'); p.green = h('#008300');
p.magenta = h('#e87ba4'); p.ink = h('#0b0b0b'); p.ink2 = h('#52514e'); p.light = h('#d8d7d2');
end
