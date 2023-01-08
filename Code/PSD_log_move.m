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