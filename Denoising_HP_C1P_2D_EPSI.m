%%
% Tensor denoising for HP C1Pyruvate 2D EPSI data
% donghyun.hong@ucsf.edu
%
clear all
close all
clc

addpath(genpath(pwd));
load('example_data/C1Pyr_EPSI_shTERT.mat');

%% Zerofilling 
zf_factor = 512;
[MetImageC_TD, MetImageC_TD_ZF, MetImageC_FD_ZF] = Zero_filling(MetImageC, zf_factor);
Data_5D_orig = permute(MetImageC_TD_ZF,[3 1 2 4]);
size_Data_5D_orig = size(Data_5D_orig);
idxarry_ppm_zf = imresize(idxarry_ppm, [1 size_Data_5D_orig(1)]);
%% TD regularization
rank = 5;
tensorRank = [rank, rank, rank, size_Data_5D_orig(4)];
Data_5D_denoise = td_regularization(Data_5D_orig,tensorRank);

Data_5D_orig_FD = fft(Data_5D_orig);
Data_5D_denoise_FD = fft(Data_5D_denoise);

specNoiseRegion = 1:100;
noise_before = std(Data_5D_orig_FD(specNoiseRegion,1,1,1));
noise_after = std(Data_5D_denoise_FD(specNoiseRegion,1,1,1));
fprintf('Noise reduction = %01d \n',noise_before/noise_after);

%% plot 2D EPSI images
% specify metabolote area
spec_idx_met = {5:100,250:350,400:500}; %Pyr (C1P,ala,Lac)

met_specAUC_before = spec2img(Data_5D_orig_FD,spec_idx_met);
met_specAUC_after = spec2img(Data_5D_denoise_FD,spec_idx_met);

plot_temporal(met_specAUC_before,met_specAUC_after);

t=3;

%%Additional part
ymin1 = min(min(min(real(Data_5D_orig(:,:,:,t)))));
ymax1 = max(max(max(real(Data_5D_orig(:,:,:,t)))));
ymin2 = min(min(min(real(Data_5D_denoise(:,:,:,t)))));
ymax2 = max(max(max(real(Data_5D_denoise(:,:,:,t)))));
ymin3 = min(min(min(abs(fft(Data_5D_orig(:,:,:,t))))));
ymax3 = max(max(max(abs(fft(Data_5D_orig(:,:,:,t))))));
ymin4 = min(min(min(abs(fft(Data_5D_denoise(:,:,:,t))))));
ymax4 = max(max(max(abs(fft(Data_5D_denoise(:,:,:,t))))));

%% before and after denoising

figure()
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

[~,~,t] = size(image);

axes('Position', [.015, .05, 1, 0.9]);
imagesc(imresize(flipud(image(:,:,t-13)'),[200 256]),'AlphaData', .2) ;hold on;
axis image off; colormap(jet);

voxel_x = 5;
voxel_y = 4;

%% heatmap
[s1 s2 s3] = size(squeeze(abs(Data_5D_denoise_FD(:,:,:,3))));
[maxval, ind] = max(reshape(F(:), s2*s3, []));
[i j] = ind2sub([s2 s3], ind);
Max_amp = [maxval' i' j'];

for i = 1:size(Data_5D_denoise_FD,4); %timeframes
    for b = 1:y; %y-axis
        for a = 1:x; %x-axis
%spectra_range_Pyr = find((idxarry_ppm > ppm_a)&(idxarry_ppm < ppm_b));
%dyn_range = i;%1:n_rep;% dynamic range for integration

inim_Pyr(a,b,i) = squeeze(sum(abs(Data_5D_denoise_FD(400:500,a,b,i))));
inim_Lac(a,b,i) = squeeze(sum(abs(Data_5D_denoise_FD(20:60,a,b,i))));
inim_Ala(a,b,i) = squeeze(sum(abs(Data_5D_denoise_FD(270:320,a,b,i))));
%inim_Pyr(:,:,i) = imresize(sum(sum(abs(yiC_zf(:,:,spectra_range_Pyr,dyn_range)),3),4),[dimension(1) dimension(3)]*resize_fac);
        end
    end
end

inim_PyrLac = inim_Lac ./ inim_Pyr;

[x_pyr, y_pyr, z_pyr] = ind2sub(size(inim_Pyr), find(inim_Pyr == max(inim_Pyr(:))));
[x_lac, y_lac, z_lac] = ind2sub(size(inim_Lac), find(inim_Lac == max(inim_Lac(:))));
[x_ala, y_ala, z_ala] = ind2sub(size(inim_Ala), find(inim_Ala == max(inim_Ala(:))));
[x_PyrLac, y_PyrLac, z_PyrLac] = ind2sub(size(inim_PyrLac), find(inim_PyrLac == max(inim_PyrLac(:))));
max_amp_Pyr = inim_Pyr(x_pyr,y_pyr,z_pyr);
max_amp_Lac = inim_Lac(x_lac,y_lac,z_lac);
max_amp_Ala = inim_Ala(x_ala,y_ala,z_ala);
max_amp_PyrLac = inim_PyrLac(x_PyrLac,y_PyrLac,z_PyrLac);

inim_HPyrLac = inim_Lac ./ max_amp_Pyr;
[x_HPyrLac, y_HPyrLac, z_HPyrLac] = ind2sub(size(inim_HPyrLac), find(inim_HPyrLac == max(inim_HPyrLac(:))));
max_amp_HPyrLac = inim_HPyrLac(x_HPyrLac,y_HPyrLac,z_HPyrLac);

% figure()
% for k = 1:size(Data_5D_denoise_FD,4)
% subplot(1,size(Data_5D_denoise_FD,4),k)
% imshow(inim_Pyr(:,:,k));
% %heatmap(inim_Pyr,[],[],[],'ColorMap', @cool, 'NaNColor', [0 0 0], 'colorbar', true);
% ylim([0 ymax4]);
% axis off;
% 
% end

[~,~,t] = size(inim_Pyr);
for time_frame = 1:t;
Full_frames{time_frame} = (flipud(image(:,:,8)'));
end

Full_Met_frames=cell(1,20);
inim_Pyr_resized = fliplr(imresize(inim_Pyr,32));
inim_Lac_resized = fliplr(imresize(inim_Lac,32));
inim_Ala_resized = fliplr(imresize(inim_Ala,32));
inim_PyrLac_resized = fliplr(imresize(inim_PyrLac,32));
inim_HPyrLac_resized = fliplr(imresize(inim_HPyrLac,32));

for time_frame = 1:t;
Full_Met_frames{time_frame} = (inim_Pyr_resized(:,:,time_frame));  
end

figure() 
%ax(1) = subplot(2,1,1)
anat_gray = montage(Full_frames,'Size', [1,t/2],'DisplayRange', []); hold on;
axis image off;
%colormap(anat_gray,'gray');
alpha(anat_gray,1);

%ax(2) = subplot(2,1,2)
met_heat = montage(Full_Met_frames,'Size', [1,t/2],'DisplayRange', []); hold on;
axis image off;
%colormap(met_heat,'jet');
alpha(met_heat,0.5);

%%
image_x = 128;
image_y = 128;
heat_x = 128;
heat_y = 128;

t1 = repmat(image(image_x-127:image_x+127,image_y-127:image_y+127,11),[1 1 20]);
imagescn_overlay(rot90(t1),[0 max(t1(:))], inim_HPyrLac_resized(heat_x-127:heat_x+127,heat_y-127:heat_y+127,:), [0 0.1],[1 10],0,0.5,'jet');