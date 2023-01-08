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