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