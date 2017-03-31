% Clear workspace to initialise
clear;

% Load input
disp('Loading data...');
% Load raw data (.csv file)
% Location of raw input file
source_raw = 'Input\EFS_S00_CO2PROFILE_R_03_2017.csv';
% Create datastore to access collection of data
ds_raw = datastore(source_raw);
% Select variable of interest
ds_raw.SelectedVariableNames = {'Date','Time','ValveNo','CO2_Avg',}; 
ds_raw.TreatAsMissing = 'NA';
% Read selected variables, save it in the workspace as a table
Raw_Data = readall(ds_raw);
% Get data from the table, change format if necessary
CO2_prof = Raw_Data.CO2_Avg;
Date_prof = Raw_Data.Date;
Time_prof = Raw_Data.Time;
ValveNo = Raw_Data.ValveNo;
Date_prof = datetime(Date_prof,'InputFormat','dd/MM/yyyy'); Date_prof.Format = 'default';
Time_prof = datetime(Time_prof,'InputFormat','HH:mm:ss');
datetime_prof = Date_prof + timeofday(Time_prof);
datetime_prof = datetime_prof + duration(00,00,01); % there is a lag of 1 sec in raw data
% clear unused variables
clearvars Raw_Data ds_raw source_raw Date_prof Time_prof; 

% Load raw data (.csv file)
% Location of raw input file
source_raw = 'Input\EddyFlux_slow_met_03_2017.csv';
% Create datastore to access collection of data
ds_raw = datastore(source_raw);
% Select variable of interest
ds_raw.SelectedVariableNames = {'Date','Time','Ta_HMP_01_Avg','Ta_HMP_155_Avg','ps_7500_Avg'}; 
ds_raw.TreatAsMissing = 'NA';
% Read selected variables, save it in the workspace as a table
Raw_Data = readall(ds_raw);
% Get data from the table, change format if necessary
Date_met = Raw_Data.Date;
Time_met = Raw_Data.Time;
ta_30 = Raw_Data.Ta_HMP_155_Avg;
ta_07 = Raw_Data.Ta_HMP_01_Avg;
P = Raw_Data.ps_7500_Avg;
Date_met = datetime(Date_met,'InputFormat','dd/MM/yyyy'); Date_met.Format = 'default';
Time_met = datetime(Time_met,'InputFormat','HH:mm:ss');
datetime_met = Date_met + timeofday(Time_met);
% clear unused variables
clearvars Raw_Data ds_raw source_raw Date_met Time_met; 
disp('Input loaded.');

% we need a continuous timestamp dataset
% Timestamps refers to end of half-hour measured.
t1 = datetime_met(1); 
t2 = datetime_met(length(datetime_met));
datetime_met_c = t1:minutes(30):t2; % t_met is a continuous timestamp vector, from 1 half-hour to last half-hour of met data
datetime_met_c = datetime_met_c';
clearvars t1 t2;
% The following lines merge DateTime and t_met, to fill timestamp gaps with
% NaN
% not needed if no timestamps gaps
% n = length(datetime_met_c);
% Ta_ct = NaN(n,2);
% RH_ct = NaN(n,2);
% P_ct = NaN(n,1);
% Ts_ct = NaN(n,4);
% SWC_ct = NaN(n,4);
% [Lia,LocB] = ismember(t_met,DateTime);
% for i = 1:n
% if LocB(i)>0
% Ta_ct(i,1) = Ta_HMP_01_Avg(LocB(i),1);
% Ta_ct(i,2) = Ta_HMP_155_Avg(LocB(i),1);
% P_ct(i,1) = ps_7500_Avg(LocB(i),1);
% end
% end

% The following lines creates a timestamp vector similar to the structure
% of raw profiler data timestamp, but with no gaps
disp('Create a continuous datetime vector...');
j=1;
datetime_prof_c(1,1) = datetime_prof(1); % start time
n = length(datetime_met_c)*8; % datetime_met_c is a continuous half-hourly timestep vector, containing timestep for met and fluxes
%progressbar
while j < n 
    %progressbar(j/n)
for i = 1:7
    datetime_prof_c(j+i,1) = datetime_prof_c(j,1) + minutes(i);
end
    j = j + 8;
    datetime_prof_c(j,1) = datetime_prof_c(j-8,1) + minutes(30);
end
datetime_prof_c(length(datetime_prof_c)) = [];
clearvars j i n;
disp('Continuous datetime vector created.');
% This script merge the continuous profiler timestamp created with the raw
% profiler timestamp, to avoid timestamp gaps. Missing timestamps line will be NaN
% data.

% Determine which elements of A are also in B as well as their corresponding locations in B.
% [Lia,Locb] = ismember(A,B)
% A = [5 3 4 2]; B = [2 4 4 4 6 8];
% Lia = 0     0     1     1
% Locb = 0     0     2     1
% The lowest index to A(3) is B(2).
% A(4) is found in B(1).

% t = no gaps timesteps
% DateTime_profiler = profiler timesteps (with gaps)
% ProfileRAW = profiler raw matrix data, without datetime column
[Lia,LocB] = ismember(datetime_prof_c,datetime_prof);
n = length(datetime_met_c)*8;
CO2_prof_c = NaN(n,1);
ValveNo_c = NaN(n,1);
progressbar
for i = 1:n
    progressbar(i/n)
    if LocB(i) > 0
        CO2_prof_c(i,:) = CO2_prof(LocB(i),:);
        ValveNo_c(i,:) = ValveNo(LocB(i),:);
    end
end
clearvars i n Lia LocB;

% Some CO2 are in wrong column, to do before removing outliers !
% Data_OK is the profiler data, with no timestamp gaps (continuous
% timestamp with NaN lines if missing timestamp)
% n = size(t_met,1)*8;
% for i = 1:n
%     if Data_OK(i,4) > 1000
%         Data_OK(i,4) = Data_OK(i,5);
%     end
% end

% Plot each height of raw CO2 for manual check
% figure;
% for i = 1:8
%     hold on;
%     use = find(ValveNo_c == i);
%     plot(datetime_prof_c(use),CO2_prof_c(use));
% end
% clearvars i use;
% Note: example: bad data between 2016/12/01 06:00:00 and 2016/12/05
% 09:00:00, all height gives the same CO2

% Outliers removal
% Decide what is an outlier based on figure ploted in script 4.
disp('Remove outliers, convert unit from ppm to umol m-3...');
n = length(datetime_met_c)*8; % or, length(datetime_prof_c), which is profiler timestamps
% Note: Profiler has 8 data per half-hour, hence length(t) = length(t_met)*8
for i = 1:n
    if CO2_prof_c(i) < 350 || CO2_prof_c(i) > 620
        CO2_prof_c(i) = NaN;
    end
end
clearvars i n;

% Convert ppm in umolco2m-3 
n = length(datetime_met_c)*8;
Height = [0.5 1 2 3.5 7 12 20 29];
Height_mobile = [0.5 1 2 3.5 7 12 14];
% Height previous to 05/10/2016 : [1 5 9 13 17 21 25 29]

% Calculation of Ta vector as linear interpolation between Ta 30m and Ta7m
% Ta(Height) = (Ta30-Ta7)/(30-4)*Height + Ta30-((Ta30-Ta7)/(30-4)*30)
Ta = NaN(n,1);
for i = 1:n/8
    for h = 1:8
        Ta(8*(i-1)+h) = (ta_30(i)-ta_07(i))/(30-4)*Height(h)+ta_30(i)-((ta_30(i)-ta_07(i))/(30-4)*30);
    end
end
clearvars i n h;

% Calculation of CO2, from ppm to umolCO2m-3, using ideal gas law, pressure
% (one per half our, measured by 7500)
% and air temperature measurements (linear interp, see above)
n = length(datetime_met_c)*8;
j=1;
CO2_umolm3 = NaN(n,1);
for i = 1:n
    if i > 8*j 
        j = j+1; 
    end
    CO2_umolm3(i) = CO2_prof_c(i)/(8.3144598*(273.15+Ta(i))/(P(j)*1000)); 
end
clearvars i n j; 

% Calculation of [CO2] at 14m as linear interpolation between 12 and 20m
% [CO2]
n = length(datetime_met_c);
CO2_14 = NaN(n,1);
for i = 1:n
    slope = (CO2_umolm3((i-1)*8+7) - CO2_umolm3((i-1)*8+6))/(20-12);
    intercept = CO2_umolm3((i-1)*8+6) - slope*12;
    CO2_14(i) = intercept + slope*14;
end
clearvars i n;

% Manual check for P_ct and Ta 
% figure; plot(datetime_met_c,P);
% figure; plot(datetime_prof_c,Ta);
disp('Storage calculation...');
% Storage calculation
n = length(datetime_met);
progressbar
Storage_layer = NaN(n*8-8,1); % -8 because need 2 half-hour for 1 storage
Storage_layer_mobile14 = NaN(n*8-8,1); % using linear int [co2] at 14m
Storage_flux = NaN(n-1,1);
Storage_flux_mobile = NaN(n-1,1);
Storage_flux_understory = NaN(n-1,1);
for i = 1:n-1;
    progressbar(i/(n-1))
    % These 8 lines calculate the change in storage flux for each dheight
    % (between 0 and 0.5m, then between 0.5m and 1m, then between 1m and 2m
    % and so on)
    Storage_layer(8*i-7,1) = (CO2_umolm3(8*i+1,1)-CO2_umolm3(8*i-7,1))/1800*Height(1);
    Storage_layer(8*i-6,1) = ((CO2_umolm3(8*i+1,1)-CO2_umolm3(8*i-7,1))/1800+(CO2_umolm3(8*i+2,1)-CO2_umolm3(8*i-6,1))/1800)*(Height(2)-Height(1))/2;
    Storage_layer(8*i-5,1) = ((CO2_umolm3(8*i+2,1)-CO2_umolm3(8*i-6,1))/1800+(CO2_umolm3(8*i+3,1)-CO2_umolm3(8*i-5,1))/1800)*(Height(3)-Height(2))/2;
    Storage_layer(8*i-4,1) = ((CO2_umolm3(8*i+3,1)-CO2_umolm3(8*i-5,1))/1800+(CO2_umolm3(8*i+4,1)-CO2_umolm3(8*i-4,1))/1800)*(Height(4)-Height(3))/2;
    Storage_layer(8*i-3,1) = ((CO2_umolm3(8*i+4,1)-CO2_umolm3(8*i-4,1))/1800+(CO2_umolm3(8*i+5,1)-CO2_umolm3(8*i-3,1))/1800)*(Height(5)-Height(4))/2;
    Storage_layer(8*i-2,1) = ((CO2_umolm3(8*i+5,1)-CO2_umolm3(8*i-3,1))/1800+(CO2_umolm3(8*i+6,1)-CO2_umolm3(8*i-2,1))/1800)*(Height(6)-Height(5))/2;
    Storage_layer(8*i-1,1) = ((CO2_umolm3(8*i+6,1)-CO2_umolm3(8*i-2,1))/1800+(CO2_umolm3(8*i+7,1)-CO2_umolm3(8*i-1,1))/1800)*(Height(7)-Height(6))/2;
    Storage_layer(8*i-0,1) = ((CO2_umolm3(8*i+7,1)-CO2_umolm3(8*i-1,1))/1800+(CO2_umolm3(8*i+8,1)-CO2_umolm3(8*i-0,1))/1800)*(Height(8)-Height(7))/2;
    Storage_layer_mobile14(i,1) = ((CO2_umolm3(8*i+6,1)-CO2_umolm3(8*i-2,1))/1800+(CO2_14(i+1,1)-CO2_14(i,1))/1800)*(Height_mobile(7)-Height_mobile(6))/2;
    % this line sum the change in storage of all layers/height
    Storage_flux(i,1) = sum(Storage_layer(8*i-7:8*i,1));
    Storage_flux_mobile(i,1) = sum(Storage_layer(8*i-7:8*i-2,1)) + Storage_layer_mobile14(i,1);
    Storage_flux_understory(i,1) = sum(Storage_layer(8*i-7:8*i-5,1));
    % example: Main tower
    % 8*1-7:8*1 = [1 2 3 4 5 6 7 8]
    % 8*2-7:8*2 = [9 10 11 12 13 14 15 16]
    % and so on...
    % Mobile
    % 8*1-7:8*1-2 = [1 2 3 4 5 6]
    % 8*2-7:8*2-2 = [9 10 11 12 13 14]
end
clearvars i n;

% See Ian's paper for formula, after Yang et al. 2007
% Basicaly, it is 
% Sum (dCO2/dt * dz)
% Between 0 and first height, we have 1 CO2 measurement per half-hour, so it is straight
% forward
% Between 2 heights after the first, we have 2 CO2 measurement per half
% hour, hence we simply average the CO2 concentration ((dCO2/dt + dCO2/dt)/2) * dz)

% Timesteps storage (end of half hour)

t1 = datetime_prof_c(9); % because timestamp refers to end of averaging period
t2 = datetime_prof_c(length(datetime_prof_c)-7);
datetime_s = t1:minutes(30):t2;
datetime_s=datetime_s';
clearvars t1 t2;
disp('Calculations done!');
% figure; plot(datetime_s,Storage_flux);

% Change in storage flux:
% When profiler is good, Sc = profilerSc, qc_Sc = 0
% When 1 or 2 inlet are bad, Sc = profilerSc(inlet working), qc_Sc = 1
% When profiler is bad, Sc = 1ptestimation, qc_Sc = 2

% example:
% For 2016, we saw manually in the plot in script 4. that profiler was down
% between 2016/12/01 06:00:00 and 2016/12/05 09:00:00
% 
% n = length(Storage_flux);
% Sc = NaN(n,1);
% qc_Sc = NaN(n,1);
% progressbar
% for i = 1:n
%     progressbar(i/n)
%     if t_s(i) >= datetime(2016,10,04,00,30,00) && t_s(i) <= datetime(2016,12,01,06,00,00) || ...
%        t_s(i) >= datetime(2016,12,05,09,00,00) 
%             Sc(i) = Storage_flux(i);
%             qc_Sc(i) = 0;
%     else
%             Sc(i) = co2_strg(i); % Loaded from EddyPro output
%             qc_Sc(i) = 2;
%     end
% end

% Plot time series Sc, and time series of raw CO2
% figure; subplot(2,1,1); plot(t_s,Sc); hold on; plot([datenum(t_s(1)) datenum(t_s(length(t_s)))],[0 0]);
% subplot(2,1,2);
% for i = 1:8
%     hold on;
%     use = find(Data_OK(:,1) == i);
%     plot(DateTime_profiler(use),Data_OK(use,4));
% end

% Note: Line under are an old version, we can come back to qc_Sc = 1
% later...
% 
% Sc = NaN(43438,1);
% qc_Sc = NaN(43438,1);
% progressbar
% for i = 1:43438
%     progressbar(i/43438)
%     if t_s(i) >= datetime(2013,11,08,15,00,00) && t_s(i) <= datetime(2014,01,06,17,30,00) || ...
%        t_s(i) >= datetime(2014,01,14,10,00,00) && t_s(i) <= datetime(2014,01,16,00,00,00) || ...
%        t_s(i) >= datetime(2014,02,27,13,00,00) && t_s(i) <= datetime(2014,04,07,10,00,00) || ...
%        t_s(i) >= datetime(2014,04,10,12,00,00) && t_s(i) <= datetime(2014,11,11,13,00,00) || ...
%        t_s(i) >= datetime(2014,11,12,15,00,00) && t_s(i) <= datetime(2014,11,30,02,00,00) || ...
%        t_s(i) >= datetime(2014,12,01,15,30,00) && t_s(i) <= datetime(2014,12,02,02,00,00) || ...
%        t_s(i) >= datetime(2015,05,05,15,00,00) && t_s(i) <= datetime(2015,09,27,07,00,00) || ...
%        t_s(i) >= datetime(2015,09,30,14,30,00) && t_s(i) <= datetime(2015,10,25,05,00,00) || ...
%        t_s(i) >= datetime(2015,12,03,14,30,00) && t_s(i) <= datetime(2015,12,07,18,00,00) || ...
%        t_s(i) >= datetime(2015,12,18,14,30,00) 
%             Sc(i) = Storage_flux(i);
%             qc_Sc(i) = 0;
%     else
%             Sc(i) = co2_strg(i);
%             qc_Sc(i) = 2;
%     end
% end