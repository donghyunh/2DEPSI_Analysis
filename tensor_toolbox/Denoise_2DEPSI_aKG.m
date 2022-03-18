%%
% v0: plot single channel denoised vs original
clear all
close all
clc

addpath(genpath(pwd));

% load 2D EPSI data
% coil = 3;
% load('example_data/b`.mat');
% 
% Data_5D = permute(csi_phased(:,:,:,:,coil,:),[4 1 3 2 5 6]);
% 
% Data_5D_orig = squeeze(Data_5D);
% size_Data_5D_orig = size(Data_5D_orig);

%%test
load('example_data_aKG/20219716_BT257_C5aKG_EPSI.mat');
%load('example_data_Pyr/20210426_C1Pyr_EPSI_shTERT.mat');


zf_factor = 512;
[MetImageC_TD, MetImageC_TD_ZF, MetImageC_FD_ZF] = Zero_filling(MetImageC, zf_factor);
%size_MetImageC = size(yiC_FID);
%Data_5D_orig = permute(yiC_FID_zf,[3 1 2 4]);
Data_5D_orig = permute(MetImageC_TD_ZF,[3 1 2 4]);
size_Data_5D_orig = size(Data_5D_orig);
idxarry_ppm_zf = imresize(idxarry_ppm, [1 size_Data_5D_orig(1)]);
%% TD regularization
 rank = 3;
 %tensorRank = [rank, rank, rank, size_Data_5D_orig(4)];
 tensorRank = [1, 4, 4, 20];
 Data_5D_denoise = td_regularization(Data_5D_orig,tensorRank);


Data_5D_orig_FD = fft(Data_5D_orig);
Data_5D_denoise_FD = fft(Data_5D_denoise);

%specNoiseRegion = 900:1000;
specNoiseRegion = 10:20;
noise_before = std(Data_5D_orig_FD(specNoiseRegion,1,1,1));
noise_after = std(Data_5D_denoise_FD(specNoiseRegion,1,1,1));
fprintf('SNR improvement = %01d \n',noise_before/noise_after);

%% plot 2D EPSI images
%spec_idx_met = {180:230,400:550,1:1000}; %aKG
%spec_idx_met = {5:15,32:42,61:87}; %aKG(2HG-Glu-aKG)
%spec_idx_met = {25:1000,600:800,55:65}; %GL_Zf_1024
spec_idx_met = {30:90,205:235,460:495}; %Pyr (C1P,ala,Lac)
%spec_idx_met = {35:45,55:65,60:70}; %Glc

met_specAUC_before = spec2img(Data_5D_orig_FD,spec_idx_met);
met_specAUC_after = spec2img(Data_5D_denoise_FD,spec_idx_met);

plot_temporal(met_specAUC_before,met_specAUC_after);

t=5;

%%Additional part
ymin1 = min(min(min(real(Data_5D_orig(:,:,:,t)))));
ymax1 = max(max(max(real(Data_5D_orig(:,:,:,t)))));
ymin2 = min(min(min(real(Data_5D_denoise(:,:,:,t)))));
ymax2 = max(max(max(real(Data_5D_denoise(:,:,:,t)))));

figure(11)
for x=1:size_Data_5D_orig(2)
    for y=1:size_Data_5D_orig(3)  
        
    subplot(size_Data_5D_orig(2),size_Data_5D_orig(3)*2+2,((size_Data_5D_orig(2)*2)+2)*(x-1)+y)
    plot(real(Data_5D_orig(:,x,y,t))); hold on;
    plot(imag(Data_5D_orig(:,x,y,t)));
    ylim([ymin1 ymax1])
    axis off;
    
    subplot(size_Data_5D_orig(2),size_Data_5D_orig(3)*2+2,((size_Data_5D_orig(2)*2)+2)*(x-1)+size_Data_5D_orig(3)+y+2)
    plot(real(Data_5D_denoise(:,x,y,t))); hold on;
    plot(imag(Data_5D_denoise(:,x,y,t))); 
    ylim([ymin2 ymax2])
    axis off;
    
    
    end
end

ymin3 = min(min(min(abs(fft(Data_5D_orig(:,:,:,t))))));
ymax3 = max(max(max(abs(fft(Data_5D_orig(:,:,:,t))))));
ymin4 = min(min(min(abs(fft(Data_5D_denoise(:,:,:,t))))));
ymax4 = max(max(max(abs(fft(Data_5D_denoise(:,:,:,t))))));

figure(12)
for x=1:size_Data_5D_orig(2)
    for y=1:size_Data_5D_orig(3)  
        
    subplot(size_Data_5D_orig(2),size_Data_5D_orig(3)*2+2,((size_Data_5D_orig(2)*2)+2)*(x-1)+y)
    plot(abs(fftshift(fft(Data_5D_orig(:,x,y,t)))));
    ylim([0 ymax3])
    axis off;
    
    subplot(size_Data_5D_orig(2),size_Data_5D_orig(3)*2+2,((size_Data_5D_orig(2)*2)+2)*(x-1)+size_Data_5D_orig(3)+y+2)
    plot(abs(fftshift(fft(Data_5D_denoise(:,x,y,t)))));
    ylim([0 ymax4])
    axis off;
    
    
    end
end

[~,~,z] = size(image);

figure(13) 
for image_frame = 1:z;
Full_frames{image_frame} = (flipud(image(:,:,image_frame)'));  
end
montage(Full_frames, 'Size', [1,z],'DisplayRange', []);
axis image off; colormap(gray);

figure(14)
for x=1:8
    for y=1:8  
    subplot(8,8,8*(x-1)+y)
    plot(abs(fftshift(fft(Data_5D_denoise(:,x,y,t)))),'linewidth',1,'color','m');
    ylim([0 ymax4])
    %set(gcf, 'Color', 'None')
    axis off;
    end
end
axes('Position', [.015, .05, 1, 0.9]);
imagesc(imresize(flipud(image(:,:,z-13)'),[200 256]),'AlphaData', .2) ;hold on;
%imagesc(imresize((RefImage),[200 256]),'AlphaData', .2) ;hold on;
axis image off; colormap(jet);

voxel_x = 5;
voxel_y = 3;

figure(15)
subplot(2,1,1)                                                                                                                                                                                                                                        
stackedplot(flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_orig(:,voxel_x,voxel_y,:)))))))),1,1);
subplot(2,1,2)                                                                                                                                                                                                                                        
stackedplot(flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_denoise(:,voxel_x,voxel_y,:)))))))),1,1);
% subplot(2,1,2)
% plot(idxarry_ppm_zf,abs(fftshift(fft(sum(squeeze(squeeze(squeeze(Data_5D_denoise(:,5,3,:)))),2))))); hold on;
% plot(idxarry_ppm_zf,real(fftshift(fft(sum(squeeze(squeeze(squeeze(Data_5D_denoise(:,5,3,:)))),2))))); hold on;
% plot(idxarry_ppm_zf,imag(fftshift(fft(sum(squeeze(squeeze(squeeze(Data_5D_denoise(:,5,3,:)))),2)))));
%set(gca, 'XDir','reverse')


rep_voxel_before = flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_orig(:,voxel_x,voxel_y,:))))))));
rep_voxel_after= flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_denoise(:,voxel_x,voxel_y,:))))))));
figure(16)
subplot(2,1,1) 
plot(rep_voxel_before(:,t));
subplot(2,1,2) 
plot(rep_voxel_after(:,t));


figure(17)
subplot(2,1,1) 
plot(sum(squeeze(rep_voxel_before(:,1:5)),2));
subplot(2,1,2) 
plot(sum(squeeze(rep_voxel_after(:,1:5)),2));

%%
figure(18)                                                                                                                                                                                                                                      
Zoomed_orig = (flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_orig(:,voxel_x,voxel_y,:)))))))));                                                                                                                                                                                                                                       
Zoomed_denoise = (flipud(abs(squeeze(squeeze(fftshift(fft(fftshift(Data_5D_denoise(:,voxel_x,voxel_y,:)))))))));
subplot(2,1,1) 
stackedplot((Zoomed_orig(665:715,:)));
subplot(2,1,2)
stackedplot((Zoomed_denoise(665:715,:)));

figure(19)
bar(sum(Zoomed_denoise(660:715,:)));

