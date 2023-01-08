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