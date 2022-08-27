
%---------------------------------------------------------------------
% Code for Coccolith carbon isotope vital effect model derived in McClelland et al. 2017, Nat comms. 
% Requires cocco_iso_concs in same folder. 
%
% Please cite as:  
% McClelland, H., Bruggeman, J., Hermoso, M. et al. 
% The origin of carbon isotope vital effects in coccolith calcite. 
% Nat Commun 8, 14511 (2017). https://doi.org/10.1038/ncomms14511
%---------------------------------------------------------------------

% Define function ISOMOD_analytical: 
% arguments:  
% DIC (DIC of external medium), CA (Carbonic Anhydrase conc. in all compartments), R (cell radius), 
% d_ex (isotopic composition of DIC in medium), Mu (growth rate), PP (PIC:POC),  ph (pH in external medium)

% The analytical solution to the coupled differential equations in
% "cocco_iso_concs.m" allows the system to be solved when it is close to
% equilibrium -- this isn't possible with the linear system approach. 

function Results = ISOMOD_analytical(  DIC, ph, d_ex, ...          % env
                            R, Mu, PP, CAall, ...       % bio
                            Fit2datB)                   % consts

N_x = 1;   % N_x (no. chloroplasts)

P_Ccell = Fit2datB(1);
T_A = Fit2datB(2);
T_B = Fit2datB(3);
E_FIX = Fit2datB(4);
FN = Fit2datB(5);


% Variable Inputs
  ext_pH    =  ph;      % external pH 
  Temp      =  15;      % Temperature (deg C) 
  Sal       =  35;      % Medium Sailinity 
  E_CAL     =   1;      % Known: Calcite precip fractionation from HCO3  (Zeebe & wolfgladrow)

  kgm3      =  1E3;
  
  
  %%-----------------------------------------------------------------------
  %%    Carbonate chemistry 
  %%-----------------------------------------------------------------------
  
  % Reaction constants: 
  R_igc  =  8.3144621;
  TempK  =  Temp +273.15;
  
  % Lueker 2000 - consitent with seacarb:
  K_1star  = 10.^((-1).*((3633.86 ./ TempK) - 61.2172 + 9.67770.*log(TempK) - 0.011555.*Sal + 0.0001152.*Sal.^2));
  K_1star  =  K_1star.*kgm3;                                                                       % .* 1e3 To put in units of mol ./ (m.^3)  rather than mol ./ kg
  K_Wstar  =  exp(148.96502 - 13847.26./TempK - 23.6521.*log(TempK) + ((118.67./TempK) - 5.977 + 1.0495.*log(TempK)).*Sal.^(1./2) - 0.01615.*Sal);
  K_Wstar  =  K_Wstar.*(kgm3.^2);                                                                    % .* 1e6 To put in units of mol.^2 ./ (m.^3).^2  rather than mol.^2 ./ kg.^2
  k_p1     =  exp(1246.98 - 6.19e4./TempK -183.*log(TempK));
  k_m1     =  (k_p1./K_1star) ;                           
  k_p4     =  4.7e7 .* exp(-23200./(R_igc.*TempK));
  k_p4     =  k_p4 ./ kgm3; 
  k_m4     =  k_p4.*(K_Wstar./K_1star);
 
  % Carbonic Anhydrase catalysed hydration./dehydration: 
  SpAp     =  2.7e7 ./kgm3;    % = Kcat ./ km from Uchikawa and with factor to convert from L to m.^3
  SpAm     =  SpAp./K_1star;
  
  % pH allowed to vary between compartments, but constant across species and samples: 
  pH_x     =  7.9;
  pH_v     =  7.1;
  pH_i     =  7.0; %  pH of compartments - measured by Anning et al. 1996.

  % CA concentrations in each compartment:
  CA_x     =  CAall;
  CA_v     =  CAall;
  CA_i     =  CAall;

  % Hydrogen ./ hydroxide ion concentrations 
  H_x      =  10.^(-pH_x).*1000;       
  OH_x     =  K_Wstar./H_x;
  H_v      =  10.^(-pH_v) .* 1000;  
  OH_v     =  K_Wstar./H_v;
  H_i      =  10.^(-pH_i) .*1000;  
  OH_i     =  K_Wstar./H_i;
 
  % Compound CO2 -> HCO3 rate constants = pH dependent.  
  k_CBi    = OH_i.*k_p4 + k_p1 + SpAp.*CA_i;
  k_CBx    = OH_x.*k_p4 + k_p1 + SpAp.*CA_x;
  k_CBv    = OH_v.*k_p4 + k_p1 + SpAp.*CA_v;
  
  % Compound HCO3 -> CO2 rate constants = pH dependent.  
  k_BCi    = k_m4 + H_i.*(k_m1 + SpAm.*CA_i);
  k_BCx    = k_m4 + H_x.*(k_m1 + SpAm.*CA_x);
  k_BCv    = k_m4 + H_v.*(k_m1 + SpAm.*CA_v);
 
  % Fractionation factors - all from zeebe
  E_cb1    =  -13;
  E_cb4    =  -11;
  E_bc1    =  -22;
  E_bc4    =  -20;
  E_cbCA   =  -1;
  E_bcCA   =  -10;
  
  % Compound CO2 -> HCO3 kinetic fractionation factors = pH dependent. 
  E_CBi = (E_cb4.*OH_i.*k_p4 + E_cb1.*k_p1 + E_cbCA.*SpAp.*CA_i)  ./  (OH_i.*k_p4 + k_p1 + SpAp.*CA_i);
  E_CBx = (E_cb4.*OH_x.*k_p4 + E_cb1.*k_p1 + E_cbCA.*SpAp.*CA_x)  ./  (OH_x.*k_p4 + k_p1 + SpAp.*CA_x);
  E_CBv = (E_cb4.*OH_v.*k_p4 + E_cb1.*k_p1 + E_cbCA.*SpAp.*CA_v)  ./  (OH_v.*k_p4 + k_p1 + SpAp.*CA_v);
  
  % Compound HCO3 -> CO2 kinetic fractionation factors = pH dependent. 
  E_BCi = (E_bc4.*k_m4 + E_bc1.*H_i.*k_m1 + E_bcCA.*H_i.*SpAm.*CA_i) ./ (k_m4 + H_i.*(k_m1 + SpAm.*CA_i));
  E_BCx = (E_bc4.*k_m4 + E_bc1.*H_x.*k_m1 + E_bcCA.*H_x.*SpAm.*CA_x) ./ (k_m4 + H_x.*(k_m1 + SpAm.*CA_x));
  E_BCv = (E_bc4.*k_m4 + E_bc1.*H_v.*k_m1 + E_bcCA.*H_v.*SpAm.*CA_v) ./ (k_m4 + H_v.*(k_m1 + SpAm.*CA_v));
  
  
  %%-----------------------------------------------------------------------
  %%    Universal cell constants
  %%----------------------------------------------------------------------- 
  
  %  Inferred from literature or independent data from this study 
  
  rho            = 20e3;    % Carbon density in mol./m
  d_n_factor     = 1   ;    % Correction factor: assume carbon only fixed during the light. 
  
  x_1D_fraction  = 1.4 ;    % ratio of chloroplast (spheroid) radius to cell radius
  x_shape_factor = 6   ;    % ratio of chloroplast (spheroid) width to height
  
  v_1D_fraction  = 0.6 ;    % ratio of coccolith vesicle (spheroid) radius to cell radius - from data
  v_shape_factor = 4   ;    % ratio of coccolith vesicle (spheroid) width to height - TEM images
  
  P_fx           =  1  ;    % Permeability of chloroplast membrane (ratio to cell membrane)
  P_fv           =  1  ;    % Permeability of Coccolith Vesicle membrane (ratio to cell membrane)


  % Compartment dimensions: 
  mu     =    Mu./(24.*60.*60);     %  G./R in seconds rather than days
  r      =    R.*1e-6;
  
  % Chloroplast 
  a_x    = r.*x_1D_fraction;      % radius of oblate spheroid
  c_x    = a_x./x_shape_factor;   % thickness of oblate spheroid 
  SA_x   = N_x.*2.*pi.*(a_x.^2 + (c_x.^2./(sin(acos(c_x./a_x))).*log((1+ sin(acos(c_x./a_x)))./(cos(acos(c_x./a_x))))));
  V_x    = (4./3).*pi.*a_x.^2.*c_x;
  
  % Coccolith Vesicle 
  a_v    = r.*v_1D_fraction;
  c_v    = a_v./v_shape_factor;
  SA_v   = 2.*pi.*(a_v.^2 + (c_v.^2./(sin(acos(c_v./a_v))).*log((1+ sin(acos(c_v./a_v)))./(cos(acos(c_v./a_v))))));
  V_v    = (4./3).*pi.*a_v.^2.*c_v;
  
  % Cytosol 
  SA_cell    = 4.*pi.*r.^2 ;          %  SA cytosol in m.^2
  Vol_cell   = (4./3).*pi.*r.^3 ;      %  Cell volume - in m.^3
  V_i        = Vol_cell - V_v - V_x ;  
  
  % Carbonate chemistry at equilibrium in external medium
  dicps  = dDICs_fast(DIC, d_ex, ext_pH, Temp, Sal);   %  need to convert from mM to M (DIC) ??
  B_e    = dicps(:, 1);                            % -  converting from M to mol ./ m.^3 ??
  C_e    = dicps(:, 3);                            % -  converting from M to mol ./ m.^3 ??
  d_Be   = dicps(:, 4);
  d_Ce   = dicps(:, 6);
  
  % Calc. and fixation rates.
  F_FIX  = rho.*Vol_cell.*mu.*d_n_factor;             % Fixation rate as a function of cell size and grown rate (and proportion of day growing) in mol C ./ cell ./ s
  F_CAL  = PP.*F_FIX;


  % Membrane permeabilities  
  UtCO2    = (F_FIX+F_CAL)./(C_e.*P_Ccell.*SA_cell + B_e.*FN.*(P_Ccell.*T_A).*SA_cell);  % utilisation at base level CO2 and HCO3 permeability 
  P_Bcell  = P_Ccell.*(T_A + T_B.*UtCO2);       %   Permeability of membrane to HCO3- as a fraction of the permeability to CO2 AND cell size.  
                                              %   DIC controls transcript abundance of HCO3- transporter genes (Bach 2013) and therefore K_B  
  P_Cx     = P_fx.*P_Ccell;                    %   membrane permeabilities 
  P_Bx     = P_fx.*P_Bcell; 
  P_Cv     = P_fv.*P_Ccell;                    %   membrane permeabilities 
  P_Bv     = P_fv.*P_Bcell; 
  


  %-------------------------------------------------------------------
  %    SOLVE MODEL :
  %-------------------------------------------------------------------
  
 iso_conc_output = cocco_iso_concs(C_e, B_e, F_FIX, F_CAL,... 
                                    P_Ccell, P_Cv, P_Cx, ...
                                    P_Bcell, P_Bv, P_Bx, FN, ...
                                    SA_cell,  SA_v, SA_x, ...
                                    V_i, V_v, V_x, ...
                                    k_BCi, k_BCv, k_BCx,...
                                    k_CBi, k_CBv, k_CBx,...
                                    d_Ce, d_Be, ...
                                    E_FIX, E_CAL,...
                                    E_CBi, E_CBv, E_CBx,...
                                    E_BCi, E_BCv, E_BCx);
                                
                                
  %    assign concentrations :
                                                 
  C_i  =  iso_conc_output(:,1);
  B_i  =  iso_conc_output(:,2);
  C_x  =  iso_conc_output(:,3);
  B_x  =  iso_conc_output(:,4);
  C_v  =  iso_conc_output(:,5);
  B_v  =  iso_conc_output(:,6);
  
  
  %    assign relevant isotopic compositions :
  
  dC_x  =  iso_conc_output(:,9);
  dB_v  =  iso_conc_output(:,12);
  
 
 
 SATcv  =  500E-6 .* (B_v./10000) ./ 3.7e-9;
 

Supply = P_Ccell.*SA_cell.*C_e + P_Bcell.*SA_cell.*B_e;
Tau = F_FIX.*(1 + PP) ./ Supply ;
TauO = F_FIX ./ Supply ;

C_org_bulk = dC_x + E_FIX;
C_calcite = dB_v + E_CAL;


% return nan when the saturation state is less than 1 or if the cellular
% concentrations of CO2 are less than 0
 NaNindex = not(SATcv > 1 & min([C_i, B_i, C_x, B_x, C_v, B_v], [], 2) > 0) ;
 Results   =  [C_x, B_v, C_org_bulk, C_calcite, Tau, TauO];
 Results(NaNindex, :)  =  repmat([nan, nan, nan, nan, nan, nan], [sum(NaNindex), 1]) ; 
 end



