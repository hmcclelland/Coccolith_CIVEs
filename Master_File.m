% Master script for model developed in Claxton et al., 2022. 
% Please cite the original paper when using any part of this code. 
%

close all
format short E

% Define the number of Monte-Carlo runs:
Number_of_MC_runs = 8;
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

%% ---------- Looking at the comparison with the independent bulk record -----------
% Load in foram data 
foram_data = readtable('data/calculated_foram_data.csv');

% Filter data for values in which foram data is present
Ordered_mean_cives = mean(d13C_dataO,1, 'omitnan')';
TARGET = Ordered_mean_cives(foram_data.index_for_matlab)';

% Need to select which foram vital effec to apply, default is 0.5permil
d13c_DIC = foram_data.calculated_d13_DIC_0_5_vital_effect;

% Prepare empty vars for scores and model outputs
CE_models = nan(Number_of_MC_runs,10);
scores = nan(Number_of_MC_runs,1);

for i = 1:Number_of_MC_runs
 model = forward_13C_absolute(Monte_CO2(i,:)', Monte_RR(i,:), DC_Rs_matrix(:,:,i), DC_Mus_matrix(:,:,i)...
    ,Fit2datB, m_DIC_initial(:,:,i), m_ph_initial(:,:,i)...
    , m_CO2_initial(:,:,i), size_frac_index, amount_calcite_data_matrix(:,:,i));
 model(isnan(model)) = nan;
 model = model(:,foram_data.index_for_matlab) + d13c_DIC';
 CE_models(i,:) = mean(model(:,:),1, 'omitnan')';
end

% Penalise runs that predict NaNs
CE_models(isnan(CE_models)) = 5;

% Calculate Scores (Target- Model)^2
scores = mean((repmat(TARGET,Number_of_MC_runs,1) - CE_models).^2,2, 'omitnan');
[Highest_run,highest_run_iter] = min(scores);

% Plot Figure 1 - parameter space that best fits d13C_DIC 
% upper plot = parameter space
% lower plot = pCO2 using henrys constants (calculated from the temperature
% from Anagnoustou et al 2020. Temp record interpolated for age points in
% this study.
henrys = readtable('data/henrys_constants.csv').Var2;
figure(1)
subplot(2,1,1)
index_for_search =  (1./scores)<1;
x = mean(Monte_RR,2, 'omitnan');
y = mean(Monte_CO2,2, 'omitnan');
z = 1./scores;
m = [x(~index_for_search) y(~index_for_search) z(~index_for_search)];
m = sortrows(m,3);
[xq,yq] = meshgrid(0.5:0.05:5,0.02:0.01:0.16);
zq = griddata(m(:,1),m(:,2),m(:,3),xq,yq);
zq(isnan(zq)) = 0;
colormap(parula)
contourf(xq,yq,zq,100,'linestyle','none')
c = colorbar();
c.Label.String = 'MSE^{-1}';
hold on
%contour(xq,yq,zq,2,'k')
set(gca, 'YScale', 'log')
xlabel('Mean R_{Calc}:R_{fix}')
ylabel('Mean output CO_2 mmol/kg')
yticks([0.025,0.05,0.1,0.15])


figure(1)
ax = subplot(2,1,2);
semilogy(ages,Monte_CO2(highest_run_iter,:)*1000./henrys','linestyle','none','marker','o','MarkerSize',8,'MarkerEdgeColor','k')
hold on
semilogy(ages(5:end),movmean(Monte_CO2(highest_run_iter,5:end),3)*1000./henrys(5:end)','black')
hold on
semilogy(ages(1:5),Monte_CO2(highest_run_iter,1:5)*1000./henrys(1:5)','color','black','linestyle','--');
xlabel('Age (Ma)')
ylabel('CO_2 (ppm)')
yticks([300,500,1000,2000,3500,6000,10000])
yticklabels({'300','500','1000','2000','3500','6000','10000'})
grid
ylim([200,10000])


%%-----------------------------------------------------------------------
%                              FIGURES
% -----------------------------------------------------------------------
Monte_CO2(Monte_CO2==0) = NaN;
Monte_RR(Monte_RR ==0) = NaN;

%% -----------------------------------------------------------------------
%  Plotting the modelled vital effects 
%  -----------------------------------------------------------------------
figure(2)
movegui('southwest');
% Set colours:
COL2 = cbrewer2('qual', 'Set1', 5,'PCHIP');
% Plot measured vital effects 
grid on
labels = ["3-5\mum Measured","5-8\mum Measured","8-10\mum Measured","10-15\mum Measured","15-20\mum Measured"];
d13C_dataO = table2array(LC_dc(:,[2 4 6 8 10]))';
d13C_data_for_plot = d13C_dataO - mean(d13C_dataO,1, 'omitnan');
d13C_dataO_sigma_for_plot = table2array(LC_dc(:,[3 5 7 9 11]))';
for k = 1:5
    hold on
    scatter(ages, d13C_data_for_plot(k,:), 75, COL2(k, :),'filled',...
              'MarkerEdgeColor','k',...
              'DisplayName',labels(k))
    errorbar(ages,d13C_data_for_plot(k,:) ,d13C_dataO_sigma_for_plot(k,:),'LineStyle','None'...
            ,'HandleVisibility','off','Color','k')
hold on
end

% Determine modelled vital effects and the 95% of realisation envelope
one = NaN(Number_of_MC_runs,37);
two = NaN(Number_of_MC_runs,37);
three = NaN(Number_of_MC_runs,37);
four = NaN(Number_of_MC_runs,37);
five = NaN(Number_of_MC_runs,37);

for i = 1:Number_of_MC_runs
SSd13C = forward_13C(Monte_CO2(i,:)', Monte_RR(i,:), DC_Rs_matrix(:,:,i), DC_Mus_matrix(:,:,i)...
    ,Fit2datB, m_DIC_initial(:,:,i), m_ph_initial(:,:,i)...
    , m_CO2_initial(:,:,i), size_frac_index, amount_calcite_data_matrix(:,:,i));
one(i,:) = SSd13C(1,:);
two(i,:) = SSd13C(2,:);
three(i,:) = SSd13C(3,:);
four(i,:) = SSd13C(4,:);
five(i,:) = SSd13C(5,:);
end

one_percentiles = prctile(one,[2.5,50,97.5],1);
two_percentiles = prctile(two,[2.5,50,97.5],1);
three_percentiles = prctile(three,[2.5,50,97.5],1);
four_percentiles = prctile(four,[2.5,50,97.5],1);
five_percentiles = prctile(five,[2.5,50,97.5],1);

% Out put the data
modelled_CIVEs = array2table([one_percentiles' two_percentiles' three_percentiles' four_percentiles' five_percentiles'],...
    'VariableNames',{'3to5-2','3to5_mean','3to5+2', '5to8-2','5to8_mean','5to8+2','8to10-2','8to10_mean','8to10+2','10to15-2','10to15_mean','10to15+2','15to20-2','15to20_mean','15to20+2'});
writetable( modelled_CIVEs,'modelled_CIVEs.csv','Delimiter',',')

% Plot modelled vital effects and their uncertainty envelope
age = ages';
x = [age,fliplr(age)];
% 3-5um size fraction
yone = [one_percentiles(1,:), fliplr(one_percentiles(3,:))];   % vector of upper & lower boundaries
patch(x,yone,COL2(1, :),'FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off') 
plot(age,one_percentiles(2,:),'color',COL2(1, :),'DisplayName',"3-5\mum Model Mean")
hold on

% 5-8um size fraction
ytwo = [two_percentiles(1,:), fliplr(two_percentiles(3,:))];   % vector of upper & lower boundaries
patch(x,ytwo,COL2(2, :),'FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off') 
plot(age,two_percentiles(2,:),'color',COL2(2, :),'DisplayName',"5-8\mum Model Mean")
hold on

% 8-10um size fraction
ythree = [three_percentiles(1,:), fliplr(three_percentiles(3,:))];   % vector of upper & lower boundaries
patch(x,ythree,COL2(3, :),'FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off') 
plot(age,three_percentiles(2,:),'color',COL2(3, :),'DisplayName',"8-10\mum Model Mean")
hold on

% 10-15um size fraction
yfour = [four_percentiles(1,:), fliplr(four_percentiles(3,:))];   % vector of upper & lower boundaries
patch(x(~isnan(yfour)),yfour(~isnan(yfour)),COL2(4, :),'FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off') 
plot(age,four_percentiles(2,:),'color',COL2(4, :),'DisplayName',"10-15\mum Model Mean")
hold on

% 15-20um size fraction
yfive = [five_percentiles(1,:), fliplr(five_percentiles(3,:))];   % vector of upper & lower boundaries
patch(x,yfive,COL2(5, :),'FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off') 
plot(age,five_percentiles(2,:),'color',COL2(5, :),'DisplayName',"15-20\mum Model Mean")
hold on

% Adding annotation arrows 
%xlim([35 56])
x = [0.95,0.95];
y = [0.75,0.9];
a = annotation('textarrow',x,y,'String','+CIVEs');

x = [0.95,0.95];
y = [0.65,0.5];
a = annotation('textarrow',x,y,'String','-CIVEs');

% Adjust plot layout
hline(0,'k-')
legend('location', 'SouthEast')
xlabel('Age (Ma)')
ylabel('\Delta\delta^{13}C_{rel. to mean}')


%% -----------------------------------------------------------------------
% % Plotting RR - Figure 1
% -----------------------------------------------------------------------
figure(3)
movegui('west');

% Chiasmolithus subplot
subplot(2,2,1)
plot(3:5, RR_mean(3:5),'color','k');
hold on 
patch([3:5 fliplr(3:5)], [RR_lower_percentile(3:5)' fliplr(RR_upper_percentile(3:5)')]...
    ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')

% Coccolithus subplot 
subplot(2,2,2)
plot(1:5, RR_mean(6:10),'color','k');
hold on
patch([1:5 fliplr(1:5)], [RR_lower_percentile(6:10)' fliplr(RR_upper_percentile(6:10)')]...
    ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')

% Discoaster subplot 
subplot(2,2,3)
plot(1:5, RR_mean(11:15),'color','k');
hold on
patch([1:5 fliplr(1:5)], [RR_lower_percentile(11:15)' fliplr(RR_upper_percentile(11:15)')]...
    ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')

% Reticulofenestra subplot 
subplot(2,2,4)
plot(1:4, RR_mean(16:19),'color','k');
hold on
patch([1:4 fliplr(1:4)], [RR_lower_percentile(16:19)' fliplr(RR_upper_percentile(16:19)')]...
    ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')

% Adjust subplot style 
legend('Mean')   
titles = [" Chiasmolithus" "Coccolithus" "Discoaster" "Reticulofenestra"];
for i = 1:4
subplot(2,2,i)
title(titles(i))
xlim([1,5])
ylim([0,2.5])
grid on
xticks([1 2 3 4 5 ])
xticklabels({'3-5\mum','5-8\mum','8-10\mum','10-15\mum','15-20\mum'})
ylabel('C:P')
xtickangle(45)
end

%% -----------------------------------------------------------------------
%  CO2 plot 
%  -----------------------------------------------------------------------
figure(4)
movegui('northwest');

%plot the points 
scatter(ages', Monte_CO2(highest_run_iter,:)'./ C_0_matrix(:,:,highest_run_iter),'red','filled','DisplayName','Best fit to DIC')
hold on 

% Plot the three point moving average over minimum of million year
% intervals
scatter(ages', Monte_CO2_mean,'black','filled','DisplayName','mean values')
plot(ages(2:3)', movmean((2:3),3), 'color','k', 'LineWidth', 1 ...
    ,'DisplayName','3pt Moving Average')
plot(ages(3:5)', Monte_CO2_mean(3:5), 'color','k','LineStyle','-.', 'LineWidth', 1 ...
    ,'HandleVisibility','off')
patch([ages(2:3)',fliplr(ages(2:3)')], [movmean(Monte_CO2_lower_percentile(2:3),3)...
   ,movmean(fliplr(Monte_CO2_upper_percentile(2:3)),3)]...
   ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')
plot(ages(5:end)', movmean(Monte_CO2_mean(5:end),3), 'color','k', 'LineWidth', 1 ...
    ,'HandleVisibility','off')
patch([ages(5:end)',fliplr(ages(5:end)')], [movmean(Monte_CO2_lower_percentile(5:end),3)...
   ,movmean(fliplr(Monte_CO2_upper_percentile(5:end)),3)]...
   ,'k','FaceAlpha',0.2,'EdgeAlpha',0,'HandleVisibility','off')
hold on

% Adjust labels 
ylabel('Relative CO_2')
xlabel('Age (Ma)')
xlim([min(ages) max(ages)])
title('Fitted CO2')
legend()
grid on

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
