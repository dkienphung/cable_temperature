function Nu = ct_nusselt_cc(Ra, Pr)
%CT_NUSSELT_CC  Churchill & Chu (1975), horizontal isothermal cylinder.
Nu = (0.60 + 0.387 * Ra.^(1/6) ./ (1 + (0.559 ./ Pr).^(9/16)).^(8/27)).^2;
end
