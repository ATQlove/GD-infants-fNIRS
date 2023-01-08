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

