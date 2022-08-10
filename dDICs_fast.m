%  Inputs:   
%  1. total DIC conc
%  2. bulk DIC d13C (here set to 0)
%  3. Medium pH 
%  4. Medium Temp

%  Outputs:  
%  1. HCO3^- concentration
%  2. CO3^2- concentration
%  3. CO2 (aq) concentration 
%  4. d13C HCO3^-
%  5. d13C CO3^2-
%  6. d13C CO2 (aq)

function results = dDICs(DICconc, dC_DIC, pH, Temp, Sal)

  %-----------------------------------------------
  %   Carbonate chemistry reaction constants: 
  %-----------------------------------------------
  
  TempK   = Temp +273.15   ;

    K_1star = exp(2.83655- 2307.1266 / TempK - 1.5529413 * log(TempK)...
            - (0.207608410 + 4.0484/TempK) * Sal^(1/2)...
            + 0.0846834 * Sal - 0.00654208 * Sal^(3/2)...
            + log(1 - 0.001005*Sal));
    K_2star = exp(-9.226508- 3351.6106 / TempK - 0.2005743 * log(TempK)...
            - (0.106901773 + 23.9722/TempK) * Sal^(1/2)...
            + 0.1130822 * Sal - 0.00846934 * Sal^(3/2)...
            + log(1 - 0.001005 * Sal));

DIC = DICconc;
eqH = 10.^(-pH);  % Equilibrium H+ conc
eqA = DIC ./ (1 + K_1star ./ eqH  + K_1star.*K_2star./eqH.^2);
eqB = DIC ./ (1 + eqH./ K_1star  + K_2star./eqH);
eqC = DIC ./ (1+ eqH / K_2star + eqH.^2 ./ (K_1star.*K_2star));
  

% Isotopes: 

% e_bg = dHCO3 - dCO2(g)
e_bg = -0.1141*Temp + 10.78;
% e_cg = dCO2 - dCO2(g)
e_cg = +0.0049*Temp - 1.31;
% e_tg = dCO3 - dCO2(g)
e_tg = -0.052*Temp + 7.22;

% e_cb = dCO2 - dHCO3
e_cb = (e_cg - e_bg)./(1+e_bg*1e-3);
% e_tb = dCO3 - dHCO3
e_tb = (e_tg - e_bg)./(1+e_bg*1e-3);


dC_HCO3 = dC_DIC - e_cb.*eqA ./ DIC - e_tb.*eqC./DIC;
dC_CO3  = e_tb.*(dC_HCO3*1e-3 +1) + dC_HCO3;
dC_CO2  = e_cb.*(dC_HCO3*1e-3 +1) + dC_HCO3;

results = [eqB, eqC, eqA, dC_HCO3, dC_CO3, dC_CO2];
end














