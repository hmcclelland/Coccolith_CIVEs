%  Master script for model developed in Claxton et al., 2022. 
%  Please cite the original paper when using any part of this code. 
% 

close all
format short E

% Define the number of Monte-Carlo runs:
Number_of_MC_runs = 4;
tic
%%  -----------------------------------------------------------------------
% All variables are set up in a 3-D matrix, with the third dimension 
% representing each individual run.---*
%  -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%  Assign functions folder to path
%  -----------------------------------------------------------------------
addpath(genpath('functions/'))

%% -----------------------------------------------------------------------
%  Load initial Data...
%  -----------------------------------------------------------------------

LC_dc = readtable('data/LC_data_fast/LC_downcore_d13C_rearranged_fast_with_error.csv');
LC_refs = readtable('data/LC_data_fast/LC_refs_fast.csv');
LC_refs = sortrows(LC_refs,[1,2]);
% size_frac_index is Matrix A from the manuscript
size_frac_index = load('data/LC_data_fast/LC_downcore_d13C_index_fast.csv');

% Defining the unique species, size fractions and ages
sps = unique(LC_refs.sp);
fracs = unique(LC_refs.frac);
ages = unique(LC_dc.Age);

%% ------------------------------------------------------------------
% Define constants (as in McClelland et al. 2017 table. 2):
%---------------------------------------------------------------------
P_Ccell   =   9.3e-04 ;%          ;% Cell membrane permeability to CO2.
T_A       =   7.8e-8 / P_Ccell  ;% Background permeability of membrane to HCO3- (as a fraction of P_Ccell)
T_B       =   5.1e-6 / P_Ccell  ;% HCO3 transporter upregulation factor (permeability of membrane to HCO3-
%#    increases linearly with background utilization)
E_FIX     =   -14.3             ;% Enzymatic fractionation by RuBisCO
FN        =    2.7              ;%
Fit2datB = [P_Ccell,T_A,T_B,E_FIX,FN];


%% -----------------------------------------------------------------------
%  Establishing parameter space to explore for measured CIVEs :  
%  -----------------------------------------------------------------------
% import the measured CIVE and assign uncertainty
% d13C_dataO is Matrix G from the manuscript
d13C_dataO = table2array(LC_dc(:,[2 4 6 8 10]))';
% determine the vital effects relative to the mean
% d13C_data_rel_mean is Matrix H from the manuscript
d13C_data_rel_mean = repmat((d13C_dataO - mean(d13C_dataO,1, 'omitnan')), [1,1,Number_of_MC_runs]);
% load in the uncertainty
d13C_dataO_sigma = repmat((table2array(LC_dc(:,[3 5 7 9 11]))'),[1,1,Number_of_MC_runs]);
% create a matrix of measured data points within associated error
d13C_data_matrix_rel_mean =  normrnd(d13C_data_rel_mean,d13C_dataO_sigma);


%% -----------------------------------------------------------------------
%  Establishing parameter space to explore for carbonate carbonate chemistry :  
%  -----------------------------------------------------------------------
% Carbonate system in terms of pH and DIC (M) (the isotope model takes these), 
% at variable CO2, constant Omega calcite.
% Calculated with seacarb R package. Need to generate different set of
% DIC & pH if you make a different assmption than const Omega: 
% carbin1 - carbin 10 reflects Omega calcite varying between 1 and 10 and
% assigned to a 1000x3x10 matrix to be choosen from each monte carlo run.

carbin = NaN(1000,2,10);                                                             
for i = 1:10 
    name = strcat('data/Carb_chem/',string(i),'_carboutoc.csv');
    carb_data_i = readtable(name);
    carbin(:,:,i) = table2array(carb_data_i(:, [2 3]));                              
end

m_DIC_initial = NaN(length(carbin),1,Number_of_MC_runs);
m_ph_initial = NaN(length(carbin),1,Number_of_MC_runs);
for i = 1:Number_of_MC_runs
    carb_selection = round(unifrnd(1,10));
    m_DIC_initial(:,:,i)  =  carbin(:,1,carb_selection).*1000;                       
    %  Convert DIC parameters to right units for isotope model. NB model uses
    %  units of mol / m3 which is = mM. 
    m_ph_initial(:,:,i)  =  carbin(:,2,carb_selection);                              
end

% CO2 concentration in medium:
m_carbsys_initial  = dDICs_fast(m_DIC_initial, 0, m_ph_initial, 15, 35);
m_CO2_initial = m_carbsys_initial(:,3,:);

%% -----------------------------------------------------------------------
%  Establishing parameter space to explore for growth rate and cell size :  
%  -----------------------------------------------------------------------
% set up cell size matrix and cell size error matrix
Rs = repmat(LC_refs.Cell_Radius, [1,1,Number_of_MC_runs]);
Rs_sigma = repmat(LC_refs.Cell_Radius_sigma, [1,1,Number_of_MC_runs]);
Rs_matrix =  normrnd(Rs,Rs_sigma); 

%determine volume and allometric growth rate with associated error 
VOLs = (4/3).*pi*(Rs_matrix).^3;
Mu_in = 10.^(-0.11.* log10(VOLs)+0.1); % from Aloisi 2015
Mu_in_matrix = normrnd(Mu_in,0.1);

% Restructure matrix to correct input dimensions 
DC_Rs_matrix = repmat(permute(Rs_matrix, [2 1 3]), [length(ages), 1, 1]);
DC_Mus_matrix = repmat(permute(Mu_in_matrix,[2 1 3]), [length(ages), 1, 1]);

%% -----------------------------------------------------------------------
%  Establishing parameter space to explore for lith weights: 
%  -----------------------------------------------------------------------

% Load in calculated lith weights with associated uncertainty
% lith_weights is Matrix C from the manuscript
lith_weights = repmat(LC_refs.Lith_Weights_pg,[1,1,Number_of_MC_runs]);
lith_weights_sigma = repmat(LC_refs.Lith_Weights_pg_sigma,[1,1,Number_of_MC_runs]);
lith_weights = normrnd(lith_weights,lith_weights_sigma);
% load in the % species distribution and determine total weight of each
% species to the each size fraction
% abund_data is Matrix B from the manuscript
abund_data = repmat(table2array(LC_dc(:,12:end)),[1,1,Number_of_MC_runs]);
Lith_size_data = repmat((permute(lith_weights, [ 2 1 3])), [length(ages), 1]);
amount_calcite_data_matrix = permute(abund_data,[2 1 3]) .* permute(Lith_size_data,[2 1 3]);
%% -----------------------------------------------------------------------
%  Establishing parameter space to explore for initial CO2:
%  -----------------------------------------------------------------------
% The initial CO2 is selected from a uniform distribution between 20 and
% 100% umols. The initial Rcalc:Rfix (RR) is selected from a uniform 
% distribution between 0.5 and 10.
C_0_matrix = NaN(length(ages),1,Number_of_MC_runs); % 
C_0_values = 10.^(unifrnd(log10(20),log10(100), [Number_of_MC_runs, 1]))* 0.001;
RR_0_matrix = NaN(20,1,Number_of_MC_runs); 
RR_0_values = 10.^(unifrnd(log10(0.5),log10(10), [Number_of_MC_runs, 1]));

for i =1:Number_of_MC_runs
    C_0_matrix(:,:,i) = C_0_values(i);       
    RR_0_matrix(:,:,i) = RR_0_values(i); %
end
 
%% -----------------------------------------------------------------------
%  Creating two matrices to save each optimised RR and CO2 to. 
%  -----------------------------------------------------------------------
Monte_RR = NaN(Number_of_MC_runs,20);
Monte_CO2 = NaN(Number_of_MC_runs,length(ages));

%% -----------------------------------------------------------------------
%  Begin Monte-Carlo optimisation loop
%  -----------------------------------------------------------------------
sum_of_opti_error_matrix = NaN(Number_of_MC_runs,1);
parfor monte_iters = 1:Number_of_MC_runs
 disp(monte_iters)
 
 %% -----------------------------------------------------------------------
 %  Loading in the data for each run
 %  -----------------------------------------------------------------------
 % Measured CIVEs relative to mean 
 % This is Matrix H from the manuscript
 d13C_data = d13C_data_matrix_rel_mean(:,:,monte_iters);
 
 % Carbonate chemistry
 m_DIC = m_DIC_initial(:,:,monte_iters);
 m_ph = m_ph_initial(:,:,monte_iters);
 m_carbsys = m_carbsys_initial(:,:,monte_iters);
 m_CO2 = m_CO2_initial(:,:,monte_iters);
 
 % Lith weight and calcite amount
 % This is Matrix B x C from the manuscript
 amount_calcite_data = amount_calcite_data_matrix(:,:,monte_iters);

 % Intial CO2 selection
 co20s = C_0_matrix(:,:,monte_iters);
 RR0s =  RR_0_matrix(:,:,monte_iters);
 
 % Cell size and growth rate
 DC_Rs = DC_Rs_matrix(:,:,monte_iters);
 DC_Mus = DC_Mus_matrix (:,:,monte_iters);
 
%% -----------------------------------------------------------------------
% Iterative optimisation approach(space walking):
%-------------------------------------------------------------------------
% Initialise looping variables to change 
 co2is = co20s;
 RRis = RR0s';
% Define number of runs for the iterative optimisation. The runs usually
% stablise after 5
 nround = 10;
% initialise iteratie CO2 and RR matrices to save values to 
CO2outs = NaN(nround, length(ages));
RRouts = NaN(nround, length(RR0s)); 
% Begin the iterative optimisation loop
for i = 1:nround 
    % First optimise for Rcalc:Rfix
     [outRR,RR_min_val] = fminsearch(@(r) minfunCO2PP(co2is, r, DC_Rs, DC_Mus, Fit2datB, m_DIC, m_ph, m_CO2, size_frac_index,...
                                        amount_calcite_data, d13C_data), RRis)
     if isnan(RR_min_val)
       break
     end
    % Reassign the out Rcalc:Rfix to RRin for next run and save first attempt
    % to temporary RR matrix 
    
     RRis = outRR;
     RRouts(i,:) = outRR;
    % Optimise for CO2 with new Rcalc:Rfix values.
     [outCO2,co2_minval] = fminsearch(@(c) minfunCO2PP(c, RRis, DC_Rs, DC_Mus, Fit2datB, m_DIC, m_ph, m_CO2, size_frac_index, ...
                                      amount_calcite_data, d13C_data), co2is);
     outCO2(outCO2 <0) = 0
    % Update CO2 in and save to temporary CO2 matrix                              
     co2is = outCO2;
     CO2outs(i,:) = outCO2';  
     sum_of_opti_error_matrix(monte_iters)= co2_minval + RR_min_val;
end
% Save final values
Monte_RR(monte_iters,:) = RRouts(end,:);
Monte_CO2(monte_iters,:) = CO2outs(end,:);

end

%% Restructuring the RR and CO2 output values to calculate percentiles:
%---------- RR  -----------
% Sort output by genus 
table2 = sortrows([LC_refs(:,1:2) array2table(Monte_RR')],2);
RR_s = table2array(table2(:,3:end));
% upper and lower bounds represent 95% of model realisations
RR_percentiles = prctile(RR_s,[2.5,50,97.5],2);
RR_lower_percentile = RR_percentiles(:,1);
RR_mean = RR_percentiles(:,2);
RR_upper_percentile = RR_percentiles(:,3);

%% ---------- CO2  -----------

CO_initial = NaN(Number_of_MC_runs,length(ages));
for i = 1:Number_of_MC_runs
    CO_initial(i,:) = permute(C_0_matrix(:,:,i), [2 1 3]);
end
Monte_CO2_relative_to_initial = Monte_CO2./CO_initial;
% upper and lower bounds represent 95% of model realisations
Monte_CO2_percentiles = prctile(Monte_CO2_relative_to_initial,[2.5,50,97.5],1);
Monte_CO2_lower_percentile = Monte_CO2_percentiles(1,:);
Monte_CO2_mean = Monte_CO2_percentiles(2,:);
Monte_CO2_upper_percentile = Monte_CO2_percentiles(3,:);


%%-----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------

% write the table for calcification to photosynthetic rate
RR_output = table(table2array(table2(:,1)),table2array(table2(:,2)),RR_lower_percentile,RR_mean,RR_upper_percentile,...
        'VariableNames',{'Size_frac','Genus','RR_lower(2sd)','RR_mean','RR_upper(2sd)'});
writetable(RR_output,'RR_Output.csv','Delimiter',',')

% Convert CO2 to umols and determine 1-SD envelope 
output_co2_t = table(ages, Monte_CO2_percentiles(1,:)', Monte_CO2_percentiles(2,:)', ...
    Monte_CO2_percentiles(3,:)','VariableNames',{'age','lower(2sd)','mean','upper(2sd)'});
writetable(output_co2_t,'CO2_Output_relative.csv','Delimiter',',')
toc
