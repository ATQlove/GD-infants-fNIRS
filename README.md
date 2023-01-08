# Gender difference in functional activity of 4-months-old infants during sleep: A Functional Near-Infrared Spectroscopy Study


## Overview

After deciding on the topic of analyzing gender differences in the dataset, we first familiarized ourselves with the format of the data collected in the dataset, and selected the most user-friendly format **.mat**. The data set was divided** into two categories: male and female. Initially, we intended to use four methods: PSD, brain connectivity, time series mean amplitude, and phase difference, and to analyze the differences from three perspectives: HbO, HbR, and HbT, with analysis of variance (ANOVA) and t-test. Finally, due to various considerations such as time and complexity, we used 3 methods** of **PSD, brain connection, and phase difference**, and started from 3 perspectives** of **HbO, HbR, and HbT**, and used **t-test** to analyze the variance. Please refer to the following points for details of each aspect.



## Data Structure

We used a file in .mat format in the dataset for analysis, which contains the following data content:

| Variable           | Description                                                  |
| ------------------ | ------------------------------------------------------------ |
| rsData.name        | Filename                                                     |
| rsData.d           | Raw intensity time series (time x channels). First half of the channels (channels 1-52) correspond to the first wavelength and the second half (channels 53-104) to the second wavelength. |
| rsData.t           | Time vector (seconds).                                       |
| rsData.s           | Information related with the experimental paradigm (not used in the current study - resting-state paradigm). |
| rsData.SD          | Information about the fNIRS optode layout configuration.     |
| rsData.aux         | Stores auxiliary measurements.                               |
| rsData.nchannels   | Total number of channels in the setup.                       |
| rsData.DPF         | Differential pathlength factor values for each wavelength.   |
| rsData.sf          | Sampling frequency (hz).                                     |
| rsData.reorder     | A vector to reorder channels from Homer2 order to the order presented in the manuscript. Channels were reordered for labelling purposes. |
| rsData.ODorigOrder | Optical density data in original order (Homer2 order).       |
| rsData.OD          | Optical density data on reordered time series.               |
| rsData.CH_remove   | Index of occipital channels in reordered data.               |
| rsData.clean_data  | Initial selection of clean segment to analyse. Excluding noisy periods at the beginning/end of the recording. |
| rsData.wav         | Data after applying the wavelet despike algorithm (reordered data). |
| rsData.conc        | Concentration (HbO, HbR and HbT)  data (reordered data).     |
| rsData.index       | Marks the beginning of the selected cleanest segment when reducing data duration to 5000 samples. |
| rsData.filt        | Filtered data (reordered data).                              |
| rsData.regressors  | Mean time series of HbO and HbR signals (occipital channels excluded, reordered data) |
| rsData.GSR_oxy     | HbO data after global signal regression. Data used for analysis (reordered data). |
| rsData.GSR_deoxy   | HbR data after global signal regression. Data used for analysis (reordered data). |
| rsData.robust_mat  | Adjacency matrix computed using a robust correlation approach (occipital channels excluded, reordered data). |



## Information about the infant whose signal was taken

Information about the acquired signal is recorded in the dataset, including the following types of information: **ID, gender, language environment, and age in days**, as described in the Appendix [csv][1].



## Code

### Code structure

+ All data processing is done in Matlab 2017b, encoded in UTF-8

  + read_data.m: responsible for reading in the data for male and female infants separately from the folders categorized as such, preprocessing them by formatting, removing occipital channels, Fourier transform, etc., and storing them in PSD_log.mat and indiv_fft.mat.
    + PSD_log.mat: stores the read and pre-processed average data.
    + indiv_fft.mat: holds the read and preprocessed individual data.
  + color_topplot.m: brain activation plotting using data from PSD_log.mat.
  + PSD_log_contrast.m: PSD curve plotting using data from PSD_log.mat.
  + PSD_log_move.m: plotting of PSD curves after translation using the data in PSD_log.mat.
  + PSD_data_channel_contrast.m: use the data in indiv_fft.mat to determine the difference of the PSD values of each channel and store the information of the channels with differences in xlsx files.
  + Phase_Contrast.m: Process the data, plot the Color Map of each adjacency matrix, and plot the radar map of phase difference.

![Matlab_script](.\File\Matlab_script.svg)

### read_data.m Data read-in

The processing of the code can be roughly divided into the following parts.

1. acquisition of the data path.
2. data read-in.
3. removal of the occipital channels.
4. doing Fourier transform for each individual's data_HbO equivalents and storing them in the indiv_fft.mat file to prepare for the subsequent calculation in PSD_log_contrast.m.
5. Average the data_HbO equivalents of male and female infants separately, do Fourier transform and store the processed data in the PSD_log.mat file.

```matlab
% What this file does is read the data and do some preprocessing of the data for subsequent chart presentation

% Gets the path to the data set
namelist_f = dir('D:\fNIRS\DataBase\osfstorage-archive\mat\female\*.mat');
namelist_M = dir('D:\fNIRS\DataBase\osfstorage-archive\mat\male\*.mat');
disp("the program is already running, please wait a minute") 

% initialize data_HbO and data_HbR

% Reading data from female infants
len_f = length(namelist_f);
dataF_HbO = zeros(5000,52,len_f);
dataF_HbR = zeros(5000,52,len_f);
dataF_HbT = zeros(5000,52,len_f);
for i = 1:len_f
    eval(['s',num2str(i),'=','struct;']);
    temp = eval(['s',num2str(i)]);
    file_name{i}=namelist_f(i).name;
    path = "D:\fNIRS\DataBase\osfstorage-archive\mat\female\"+file_name{i};
    temp= load(path);       %read every .mat file of female 
    temp_HbO = permute(temp.rsData.conc(:,1,:) , [1,3,2]);
    temp_HbR = permute(temp.rsData.conc(:,2,:) , [1,3,2]);
    temp_HbT = permute(temp.rsData.conc(:,3,:) , [1,3,2]);
    dataF_HbO(:,:,i) = temp_HbO;
    dataF_HbR(:,:,i) = temp_HbR;
    dataF_HbT(:,:,i) = temp_HbT;
end

%Remove data from occipital channels
dataF_HbO(:,24:26,:) = []; %channels 24:26
dataF_HbO(:,47:49,:) = []; %channels 50:52
dataF_HbR(:,24:26,:) = [];
dataF_HbR(:,47:49,:) = []; 
dataF_HbT(:,24:26,:) = [];
dataF_HbT(:,47:49,:) = []; 

% Reading data from male infants
len_m = length(namelist_M);
dataM_HbO = zeros(5000,52,len_m);
dataM_HbR = zeros(5000,52,len_m);
dataM_HbT = zeros(5000,52,len_m);
for i = 1:len_m
    eval(['s',num2str(i),'=','struct;']);
    temp = eval(['s',num2str(i)]);
    file_name{i}=namelist_M(i).name;
    path = "D:\fNIRS\DataBase\osfstorage-archive\mat\male\"+file_name{i};
    temp= load(path);       %read every .mat file of female 
    temp_HbO = permute(temp.rsData.conc(:,1,:) , [1,3,2]);
    temp_HbR = permute(temp.rsData.conc(:,2,:) , [1,3,2]);
    temp_HbT = permute(temp.rsData.conc(:,3,:) , [1,3,2]);
    dataM_HbO(:,:,i) = temp_HbO;
    dataM_HbR(:,:,i) = temp_HbR;
    dataM_HbT(:,:,i) = temp_HbT;
end

%Remove data from occipital channels
dataM_HbO(:,24:26,:) = []; %channels 24:26
dataM_HbO(:,47:49,:) = []; %channels 50:52
dataM_HbR(:,24:26,:) = [];
dataM_HbR(:,47:49,:) = []; 
dataM_HbT(:,24:26,:) = [];
dataM_HbT(:,47:49,:) = []; 

N = 5000; %Number of samples (i.e., length of data set)
sf = 8.9300; %sampling frequency (8.93Hz in this example)

% ------------------------------------------------------
% Prepare for the PSD_data_channel_contrast.m file. 
% Do Fourier transform on data_HbO for each one.
% ------------------------------------------------------
indiv_fft_HbO_f = zeros(2501,46,len_f);
indiv_fft_HbR_f = zeros(2501,46,len_f);
indiv_fft_HbT_f = zeros(2501,46,len_f);

for nf = 1:len_f
    % Compute Fourier Transform of HbO and HbR
    temp_HbO_f = fft(dataF_HbO(:,:,nf)); % Fourier Transform HbO data
    indiv_fft_HbO_f(:,:,nf) = 2*abs(temp_HbO_f(1:N/2+1, :)); % Keep only first half
    temp_HbR_f = fft(dataF_HbR(:,:,nf)); % Fourier Transform HbR data
    indiv_fft_HbR_f(:,:,nf) = 2*abs(temp_HbR_f(1:N/2+1, :)); % Keep only first half
    temp_HbT_f = fft(dataF_HbT(:,:,nf)); % Fourier Transform HbT data
    indiv_fft_HbT_f(:,:,nf) = 2*abs(temp_HbT_f(1:N/2+1, :)); % Keep only first half
end


indiv_fft_HbO_m = zeros(2501,46,len_m);
indiv_fft_HbR_m = zeros(2501,46,len_m);
indiv_fft_HbT_m = zeros(2501,46,len_m);

for nm = 1:len_m
% Compute Fourier Transform of HbO and HbR
    temp_HbO_m = fft(dataM_HbO(:,:,nm)); % Fourier Transform HbO data
    indiv_fft_HbO_m(:,:,nm) = 2*abs(temp_HbO_m(1:N/2+1, :)); % Keep only first half
    temp_HbR_m = fft(dataM_HbR(:,:,nm)); % Fourier Transform HbR data
    indiv_fft_HbR_m(:,:,nm) = 2*abs(temp_HbR_m(1:N/2+1, :)); % Keep only first half
    temp_HbT_m = fft(dataM_HbT(:,:,nm)); % Fourier Transform HbT data
    indiv_fft_HbT_m(:,:,nm) = 2*abs(temp_HbT_m(1:N/2+1, :)); % Keep only first half
end

save("indiv_fft.mat","indiv_fft_HbO_f","indiv_fft_HbR_f","indiv_fft_HbT_f","indiv_fft_HbO_m","indiv_fft_HbR_m","indiv_fft_HbT_m");


% -----------------------------------------------------------------------
% Fourier transform of data_HbO averages for male and female infants.
% -----------------------------------------------------------------------

% Find the average value of HbO, HbR, HbT data
dataF_HbO = mean(dataF_HbO,3);
dataF_HbR = mean(dataF_HbR,3);
dataF_HbT = mean(dataF_HbT,3);
dataM_HbO = mean(dataM_HbO,3);
dataM_HbR = mean(dataM_HbR,3);
dataM_HbT = mean(dataM_HbT,3);

% Calculate frequency range
freq = linspace(0, sf/2, N/2+1);
% Compute Fourier Transform of HbO and HbR
fft_HbO_f = fft(dataF_HbO); % Fourier Transform HbO data
fft_HbO_f = 2*abs(fft_HbO_f(1:N/2+1, :)); % Keep only first half
fft_HbR_f = fft(dataF_HbR); % Fourier Transform HbR data
fft_HbR_f = 2*abs(fft_HbR_f(1:N/2+1, :)); % Keep only first half
fft_HbT_f = fft(dataF_HbT); % Fourier Transform HbT data
fft_HbT_f = 2*abs(fft_HbT_f(1:N/2+1, :)); % Keep only first half

fft_HbO_m = fft(dataM_HbO); % Fourier Transform HbO data
fft_HbO_m = 2*abs(fft_HbO_m(1:N/2+1, :)); % Keep only first half
fft_HbR_m = fft(dataM_HbR); % Fourier Transform HbR data
fft_HbR_m = 2*abs(fft_HbR_m(1:N/2+1, :)); % Keep only first half
fft_HbT_m = fft(dataM_HbT); % Fourier Transform HbT data
fft_HbT_m = 2*abs(fft_HbT_m(1:N/2+1, :)); % Keep only first half

% Store the data for later use
save("PSD_log.mat","fft_HbO_f","fft_HbR_f","fft_HbT_f","fft_HbO_m","fft_HbR_m","fft_HbT_m");
```



### PSD_log_contrast.m PSD curve plotting

This code block is intended to use to draw PSD curves using the read-in data and is structured as follows.

1. setting each parameter of the generated image.
2. smoothing the curve to be plotted in order to make the variance easier to display, here the sliding average smoothing is chosen and the window size is chosen to be 50.
3. plotting the curve of the PSD, using a linear scale for the horizontal axis and a logarithmic scale for the vertical axis for ease of presentation.

```matlab
% Create figure
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', [1 1 1]);
set(0, 'DefaultAxesFontSize', 24, 'DefaultAxesTitleFontWeight', 'normal') 

% Read in data from PSD_log.mat
load("PSD_log.mat");

N = 5000; %Number of samples (i.e., length of data set)
sf = 8.9300; %sampling frequency (8.93Hz in this example)

% Calculate frequency range
freq = linspace(0, sf/2, N/2+1);

% ----------------------------------------HbO block----------------------------------------------------
% HbO female
% --------------- Visualization method (mean PSD across channels ,log scale)
smth_fft_HbO_f = smooth(mean(fft_HbO_f,2),50);
plot(freq(100:2501), smth_fft_HbO_f(100:2501), 'b', 'linewidth', 2.0); box off
set(gca, 'XScale', 'linear'); set(gca, 'YScale', 'log')
xlim([0 sf/2]); xlabel('Frequency (Hertz)')
title ('Mean PSD HbO (log scale)')
hold on

% --------------- 
% HbO male
% --------------- 
smth_fft_HbO_m = smooth(mean(fft_HbO_m,2),50);
plot(freq(100:2501), smth_fft_HbO_m(100:2501), 'r', 'linewidth', 2.0); box off

% ----------------------------------------HbR block----------------------------------------------------
% HbR female
% % --------------- Visualization method (mean PSD across channels ,log scale)
smth_fft_HbR_f = smooth(mean(fft_HbR_f,2),50);
plot(freq(200:2501), smth_fft_HbR_f(200:2501), 'color',[0.00 0.79 0.34], 'linewidth', 2.0); box off
hold on

% --------------- 
% HbR male
% --------------- 
smth_fft_HbR_m = smooth(mean(fft_HbR_m,2),50);
plot(freq(200:2501), smth_fft_HbR_m(200:2501), 'color',[0.63 0.13 0.94], 'linewidth', 2.0); box off
set(gca, 'XScale', 'linear'); set(gca, 'YScale', 'log')
hold off
```



### PSD_log_move.m PSD curve plotting for panning

Since, after plotting the PSD graph, it was observed that the PSD level values were almost identical between men and women for a certain frequency difference before and after, the following post-panning curves were plotted, driven by this initial observation, along with the output of the magnitude of the values of this difference:

```matlab
% Create figure
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', [1 1 1]);
set(0, 'DefaultAxesFontSize', 24, 'DefaultAxesTitleFontWeight', 'normal') 

% Read in data from PSD_log.mat
load("PSD_log.mat");

N = 5000; %Number of samples (i.e., length of data set)
sf = 8.9300; %sampling frequency (8.93Hz in this example)

% Calculate frequency range
freq = linspace(0, sf/2, N/2+1);

% ----------------------------------------HbO block----------------------------------------------------
% HbO female
% --------------- Visualization method (mean PSD across channels ,log scale)
smth_fft_HbO_f = smooth(mean(fft_HbO_f,2),50);
plot(freq(100:2501), smth_fft_HbO_f(100:2501), 'b', 'linewidth', 2.0); box off
set(gca, 'XScale', 'linear'); set(gca, 'YScale', 'log')
xlim([0 sf/2]); xlabel('Frequency (Hertz)')
title ('Mean PSD HbO (log scale)')
hold on

% +++++++++++++++
% HbO male 
% move left
% +++++++++++++++
x = freq(100:2501)-0.05;
plot(x, smth_fft_HbO_m(100:2501), 'g', 'linewidth', 2.0); box off
set(gca, 'XScale', 'linear'); set(gca, 'YScale', 'log')

legend('female','male')
hold off

% find out  the max and min in the curve
min_f = min(mean(fft_HbO_f,2));
min_m = min(mean(fft_HbO_m,2));
max_f = max(mean(fft_HbO_f,2));
max_m = max(mean(fft_HbO_m,2));
min_diff = min_f-min_m
max_diff = max_f-max_m
```





### PSD_data_channel_contrast.m t-test between channels

This code block performs a t-test on Fourier-transformed data for everyone.

1. we first discard the first 10 collection points (2501 collection points in total) due to collection noise.
2. since there were 51 female infant data in the dataset and only 48 male infant data, we randomly selected 48 female infant data for the alignment of the t-test (since the numbers in the dataset were not sorted according to any order, we directly randomly took 3:50 of the data here).
3. a t-test was performed on the data that underwent the above treatment, and those that passed the t-test were shown with their pass results, as well as p-values.
4. finally save all values in an xlsx file for visualization later.

```matlab
% Read in data from indiv_fft.mat
load('indiv_fft.mat');

ch = 46;  % Number of channels
output_f = [];
output_m = [];
output = [];

for nch = 1:ch
%     Discard the data from the first 10 sampling points 
%     (because this value is too large compared to the others and will affect the average), 
%     take all channels in indiv_fft_HbR_f, and 48 random female infants
        cont_f = indiv_fft_HbR_f(10:2501,nch,3:50);  
        cont_f = squeeze(cont_f);

        cont_m = indiv_fft_HbR_m(10:2501,nch,:);
        cont_m = squeeze(cont_m);

        mean_cont_f = mean(cont_f,1);
        mean_cont_f = squeeze(mean_cont_f);
        mean_cont_m = mean(cont_m,1);
        mean_cont_m = squeeze(mean_cont_m);
        
        % Perform a t-test
        [h,p] = ttest2(mean_cont_f,mean_cont_m);
        
        if(h == 1)
            nch
            h,p
            disp("----------------------------------")
        end
        
        output_f = [output_f,mean_cont_f'];
        output_m = [output_m,mean_cont_m'];
        temp = [mean_cont_f',mean_cont_m'];
        output = [output,temp];
        
end

xlswrite('parameter_f.xlsx',output_f);% Store the data for later use
xlswrite('parameter_m.xlsx',output_m);% Store the data for later use
% xlswrite('mean_cont_all.xlsx',output);
```



### Phase_Contrast.m Adjacency matrix and phase difference plotting

This code block carries out the drawing of the Color Map of the adjacency matrix and the radar map of the phase difference, and the final image is formatted and adjusted by software such as visio.

The overall steps are as follows.

1. read in the data (most of the body of the read-in data is the same as read_data.m. The only difference is that the variable rsData.filt is used in the dataset, because it is only by using the preliminary filtered data in the dataset that the values between 90⁰~270⁰ can be filtered, in line with the HbO-HbR phase difference (definition).
2. use the corr function to obtain the adjacency matrix of HbO-HbR and plot the Color Map.
3. repeating the above process.
4. find the value of hPod using the Hebelt transform, etc., and plot the radar map.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filt version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The part of the data read here is roughly the same as before, the only change is that the data read is rsData.filt in the source data set. 
% This is done to show better phase difference results by using the data that has been preliminarily filtered in the source data set.

namelist_f = dir('D:\fNIRS\DataBase\osfstorage-archive\mat\female\*.mat');
namelist_M = dir('D:\fNIRS\DataBase\osfstorage-archive\mat\male\*.mat');
% namelist format:
% name -- filename
% date -- modification date
% bytes -- number of bytes allocated to the file
% isdir -- 1 if name is a directory and 0 if not

%initialize data_HbO and data_HbR
disp("the program is already running, please wait a minute") 

len_f = length(namelist_f);
dataF_HbO = zeros(5000,52);
dataF_HbR = zeros(5000,52);
dataF_HbT = zeros(5000,52);
for i = 1:len_f
    eval(['s',num2str(i),'=','struct;']);
    temp = eval(['s',num2str(i)]);
    file_name{i}=namelist_f(i).name;
    path = "D:\fNIRS\DataBase\osfstorage-archive\mat\female\"+file_name{i};
    temp= load(path);       %read every .mat file of female 
    temp_HbO = temp.rsData.filt(:,1:52);
    temp_HbR = temp.rsData.filt(:,53:104);
%     temp_HbO = temp.rsData.GSR_oxy(:,:);
%     temp_HbR = temp.rsData.GSR_deoxy(:,:);
    temp_HbT = temp_HbR+temp_HbO;
    dataF_HbO = dataF_HbO+temp_HbO;
    dataF_HbR = dataF_HbR+temp_HbR;
    dataF_HbT = dataF_HbT+temp_HbT;
end

dataF_HbO(:,24:26) = [];
dataF_HbO(:,47:49) = []; %dataF_HbO(:,50:52) = [];
dataF_HbR(:,24:26) = [];
dataF_HbR(:,47:49) = []; %dataF_HbR(:,50:52) = [];
dataF_HbT(:,24:26) = [];
dataF_HbT(:,47:49) = []; %dataF_HbT(:,50:52) = [];


len_m = length(namelist_M);
dataM_HbO = zeros(5000,52);
dataM_HbR = zeros(5000,52);
dataM_HbT = zeros(5000,52);
for i = 1:len_m
    eval(['s',num2str(i),'=','struct;']);
    temp = eval(['s',num2str(i)]);
    file_name{i}=namelist_M(i).name;
    path = "D:\fNIRS\DataBase\osfstorage-archive\mat\male\"+file_name{i};
    temp= load(path);       %read every .mat file of female 
    temp_HbO = temp.rsData.filt(:,1:52);
    temp_HbR = temp.rsData.filt(:,53:104);
%     temp_HbO = temp.rsData.GSR_oxy(:,:);
%     temp_HbR = temp.rsData.GSR_deoxy(:,:);
    temp_HbT = temp_HbR+temp_HbO;
    dataM_HbO = dataM_HbO+temp_HbO;
    dataM_HbR = dataM_HbR+temp_HbR;
    dataM_HbT = dataM_HbT+temp_HbT;
end

dataM_HbO(:,24:26) = []; %channels 24:26
dataM_HbO(:,47:49) = []; %channels 50:52
dataM_HbR(:,24:26) = [];
dataM_HbR(:,47:49) = []; 
dataM_HbT(:,24:26) = [];
dataM_HbT(:,47:49) = []; 


dataF_HbO = dataF_HbO/len_f;
dataF_HbR = dataF_HbR/len_f;
dataF_HbT = dataF_HbT/len_f;
dataM_HbO = dataM_HbO/len_m;
dataM_HbR = dataM_HbR/len_m;
dataM_HbT = dataM_HbT/len_m;
ch = 46;


%-------------------------------Color map of adjacency matrix--------------------------------

% Create figure
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', [1 1 1]);
set(0, 'DefaultAxesFontSize', 24, 'DefaultAxesTitleFontWeight', 'normal') 

% Compute and plot adjacency matrices
% -------------------HbO---------------------
subplot(1,2,1)
imagesc(corr(dataF_HbO), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbO female'); axis square

subplot(1,2,2)
imagesc(corr(dataM_HbO), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbO male'); axis square
colormap jet

% -------------------HbR---------------------
figure(2)
subplot(1,2,1)
imagesc(corr(dataF_HbR), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbR female'); axis square

subplot(1,2,2)
imagesc(corr(dataM_HbR), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbR male'); axis square
colormap jet

% -------------------HbT---------------------
figure(3)
subplot(1,2,1)
imagesc(corr(dataF_HbT), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbT female'); axis square

subplot(1,2,2)
imagesc(corr(dataM_HbT), [-1 1]); colorbar
xlabel('Channels'); ylabel('Channels')
title('HbT male'); axis square
colormap jet

% -------------------HbO-HbR---------------------
figure(4)
subplot(1,2,1)
imagesc(corr(dataF_HbO, dataF_HbR), [-1 1]); colorbar
xlabel('Channels');
title('HbO-HbR female'); axis square

subplot(1,2,2)
imagesc(corr(dataM_HbO, dataM_HbR), [-1 1]); colorbar
xlabel('Channels');
title('HbO-HbR male'); axis square
colormap jet


%-------------------------------Radar map of phase difference--------------------------------
%Compute hPod_f value (HbO-HbR phase difference as described in Watanabe et al., 2017)
hPod_f = zeros(1, ch);
for nch = 1:ch
    
    % Calculate Hilbert transformation of the signals
    HbO_hilbert = hilbert(dataF_HbO(:, nch));
    HbR_hilbert = hilbert(dataF_HbR(:, nch));
    
    % Calculate instantaneous phase
    HbO_inst = unwrap(angle(HbO_hilbert));
    HbR_inst = unwrap(angle(HbR_hilbert));
    
    % Calculate phase difference
    ph_dif = HbO_inst - HbR_inst;
    
    % Compute and store hPod_f
    hPod_f (nch) = angle(mean(exp(sqrt(-1)*ph_dif)));
end

% Plot hpod
figure(5)
subplot(1,2,1)
h = polarhistogram (hPod_f, 10); % adjust number of bins
set(h, 'linewidth', 1, 'FaceColor', 'b')
title('Phase difference HbO-HbR (degrees)')

hPod_m = zeros(1, ch);
for nch = 1:ch
    
    % Calculate Hilbert transformation of the signals
    HbO_hilbert = hilbert(dataM_HbO(:, nch));
    HbR_hilbert = hilbert(dataM_HbR(:, nch));
    
    % Calculate instantaneous phase
    HbO_inst = unwrap(angle(HbO_hilbert));
    HbR_inst = unwrap(angle(HbR_hilbert));
    
    % Calculate phase difference
    ph_dif = HbO_inst - HbR_inst;
    
    % Compute and store hPod_f
    hPod_m (nch) = angle(mean(exp(sqrt(-1)*ph_dif)));
end

% Plot hpod
subplot(1,2,2)
h = polarhistogram (hPod_m, 10); % adjust number of bins
set(h, 'linewidth', 1, 'FaceColor', 'b')
title('Phase difference HbO-HbR (degrees)')
```



### color_topplot.m Plot brain activation map

This block of code is intended to plot the brain activation plot using the data previously read in, with the final image formatted and adjusted by software such as visio.

```matlab
% Read in data from PSD_log.mat
load("PSD_log.mat");

%S et two 1*52 0 matrices 52 refers to the number of channels
ch_x=zeros(1,52);
ch_y=zeros(1,52);

ch_x(1:10)=2:2:20;
ch_x(11:21)=1:2:21;
ch_x(22:31)=2:2:20;
ch_x(32:42)=1:2:21;
ch_x(43:52)=2:2:20;
ch_y(1:10)=5;
ch_y(11:21)=4;
ch_y(22:31)=3;
ch_y(32:42)=2;
ch_y(43:52)=1;

ch_select_p1 = 1:23;
ch_select_p2 = 27:49;

[xq,yq] = meshgrid(0:0.1:22,0:0.1:6);    %Laying 3D grid

chx_p1 = ch_x(ch_select_p1);  chy_p1 = ch_y(ch_select_p1);
chx_p2 = ch_x(ch_select_p2);  chy_p2 = ch_y(ch_select_p2);
chx_p = [chx_p1,chx_p2];   chy_p = [chy_p1,chy_p2];

showF_HbO = mean(dataF_HbO,1);
showM_HbO = mean(dataM_HbO,1);
showF_HbR = mean(dataF_HbR,1);
showM_HbR = mean(dataM_HbR,1);
showF_HbT = mean(dataF_HbT,1);
showM_HbT = mean(dataM_HbT,1);

valq1 = griddata(chx_p,chy_p,showF_HbO,xq,yq,'v4');    %Interpolation of scattered data
valq2 = griddata(chx_p,chy_p,showM_HbO,xq,yq,'v4');
valq3 = griddata(chx_p,chy_p,showF_HbR,xq,yq,'v4');    %Interpolation of scattered data
valq4 = griddata(chx_p,chy_p,showM_HbR,xq,yq,'v4');
valq5 = griddata(chx_p,chy_p,showF_HbT,xq,yq,'v4');    %Interpolation of scattered data
valq6 = griddata(chx_p,chy_p,showM_HbT,xq,yq,'v4');


%subplot(3,2,6)
pcolor(xq,yq,valq6);
shading flat;
axis off
caxis([-0.5,1.5]);
title('Tot-Hb male')
colorbar('southoutside')
set(gca,'FontSize',15)

colorbar('southoutside')
colormap('jet')
```



## Statistical analysis and its visualization

Statistical analysis and its visualization are performed using GraphPad Prism 9 software as follows.

1. first create the Project, here using Mean, SD, N benchmarks.

![temp1](https://github.com/ATQlove/Gender-difference-in-functional-activity-of-4-months-old-infants-during-sleep-A-fNIRS-study/raw/master/File/temp1.png)

1. Create multiple forms and pick out the channels from the PSD_data_channel_contrast.m file that yield discrepancies from the parameter_f.xlsx and parameter_m.xlsx files and fill in the forms at:

   ![temp2](https://github.com/ATQlove/Gender-difference-in-functional-activity-of-4-months-old-infants-during-sleep-A-fNIRS-study/raw/master/File/temp2.jpg)

2. Click New Analysis to create an analysis and use Unpaired t-test:

   ![temp3](https://github.com/ATQlove/Gender-difference-in-functional-activity-of-4-months-old-infants-during-sleep-A-fNIRS-study/raw/master/File/temp3.jpg)

3. Visualization of results plotted against journal requirements.

   ![temp4](https://github.com/ATQlove/Gender-difference-in-functional-activity-of-4-months-old-infants-during-sleep-A-fNIRS-study/raw/master/File/temp4.jpg)

4. Create layout and arrange format to meet SCI requirements.

   ![temp5](https://github.com/ATQlove/Gender-difference-in-functional-activity-of-4-months-old-infants-during-sleep-A-fNIRS-study/raw/master/File/temp5.jpg)





There is another key diagram in visualization, that is, the placement of the channel when collecting data. Find the location diagram of the 10-20 System used in the data set, combined with Visio and other software to draw.



The following data set was used for this experiment, and we sincerely thank all those who helped with this study:

Blanco B, Molnar M, Carreiras M, Caballero-Gaudes C. Rs_4months (2022) [updated 2022]. Available from: https://github.com/borjablanco/RS_4months.



[1]:./BCBL_RS4_participant-info.csv
