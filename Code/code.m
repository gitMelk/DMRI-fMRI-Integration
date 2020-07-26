%%
clc
clear all 
close all
%% Data preproccessing 
% Load data
fmri_raw = load_untouch_nii('./Data/FMRI/FMRI_2_T1_2mm.nii')
fmri = double(fmri_raw.img);
gray_matter_raw = load_untouch_nii('./Data/Structural/c1T1_2mm.nii')
white_matter_raw = load_untouch_nii('./Data/Structural/c2T1_2mm.nii')
cspinal_fluid_raw = load_untouch_nii('./Data/Structural/c3T1_2mm.nii')
GM = double(gray_matter_raw.img) ./ 255;
WM = double(white_matter_raw.img) ./ 255;
CSF = double(cspinal_fluid_raw.img) ./ 255;
% thresholding 
thresh_GM = 0.9;
thresh_WM = 0.9;
thresh_CSF = 0.8;
GM_tr = GM > thresh_GM; % threshold probaility 
WM_tr = WM > thresh_WM;
CSF_tr = CSF > thresh_CSF;
% erosion
se = [1; 1]; % cusom object, as the default one are too agressive 
GM_erode = imerode(GM_tr, se);
WM_erode = imerode(WM_tr, se);
CSF_erode = imerode(CSF_tr, se);
% Plot of original image and filtering steps
slice = 62;

subplot(3,3,1), imagesc(imrotate(squeeze(GM(:,slice,:)),90)), title('Original GM')
subplot(3,3,2), imagesc(imrotate(squeeze(WM(:,slice,:)),90)), title('Original WM')
subplot(3,3,3), imagesc(imrotate(squeeze(CSF(:,slice,:)),90)), title('Original CSF')
subplot(3,3,4), imagesc(imrotate(squeeze(GM_tr(:,slice,:)),90)), title('GM thresholded')
subplot(3,3,5), imagesc(imrotate(squeeze(WM_tr(:,slice,:)), 90)), title('WM thresholded')
subplot(3,3,6), imagesc(imrotate(squeeze(CSF_tr(:,slice,:)),90)), title('CSF thresholded')
subplot(3,3,7), imagesc(imrotate(squeeze(GM_erode(:,slice,:)),90)), title('GM thresholded')
subplot(3,3,8), imagesc(imrotate(squeeze(WM_erode(:,slice,:)),90)), title('MASK - WM eroded')
subplot(3,3,9), imagesc(imrotate(squeeze(CSF_erode(:,slice,:)),90)), title('MASK - CSF eroded')

%% Mean fMRI signal of WM and CSF

% Extraction with masks
for p=1:size(fmri,3)
    for q=1:size(fmri,4)
        WM_masked(:,:,p,q) = fmri(:,:,p,q).*WM_erode(:,:,p);
        CSF_masked(:,:,p,q) = fmri(:,:,p,q).*CSF_erode(:,:,p);
    end
end

% Mean computation
for j=1:size(fmri,4)
    tmp1 = WM_masked(:,:,:,j);
    WM_masked_mean_225(j) = mean2(tmp1);

    tmp2 = CSF_masked(:,:,:,j);
    CSF_masked_mean_225(j) = mean2(tmp2);    
end

% Plot of the results
figure
subplot(2,1,1), plot(WM_masked_mean_225), title('WM masked mean'), xlim([0 225])
subplot(2,1,2), plot(CSF_masked_mean_225), title('CSF masked mean'), xlim([0 225])

%% 1.c sumEPI
sumEPI = squeeze(sum(fmri,4));
sumEPI = squeeze(sum(sumEPI,3));

% looking at hist and the plot I decided the threshold
figure, 
subplot(221), imagesc(imrotate(sumEPI,90)), title('EPI')
sumEPI_mask = sumEPI> 0.75e7;
subplot(222), imagesc(imrotate(sumEPI.*sumEPI_mask,90)), title('MASK EPI')

CC = bwconncomp(sumEPI_mask);  % connected components
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels); % selection of the largest connected component
mask_EPI = zeros(size(sumEPI));
mask_EPI(CC.PixelIdxList{idx}) = 1;
% plot of the selected section and final mask 
subplot(223), imagesc(imrotate(sumEPI-sumEPI.*mask_EPI,90)), title('Estracted region')
subplot(224), imagesc(imrotate(sumEPI.*mask_EPI,90)), title('final MASK EPI')


%% 1.d Atlas masked with GM and sumEPI
atlas_raw = load_untouch_nii('./Data/Atlas/Hammers_2_T1_2mm_int.nii')
atlas = double(atlas_raw.img);
atlas_GM = atlas.*GM_erode;
atlas_sumEPI = atlas_GM.*mask_EPI;
tmp = 51;
figure, 
imagesc(imrotate(atlas_sumEPI(:,:,tmp),90)), title('Slice 51 - Atlas with Masks')
%% 1.e
ROI_list_discard = [3,4,17,18,19,44,45,46,47,48,49,74,75];
% ROI 43 is absent 

% remove ROIS
for j=1:size(atlas_sumEPI,3)
    temp_slice = squeeze(atlas_sumEPI(:,:,j));
    temp_mask = ones(size(temp_slice)); 
    for i=1:max(temp_slice(:))
        if (find(ROI_list_discard==i)>0)
            temp_mask(find(temp_slice==i)) = 0;            
        
        elseif(size(find(temp_slice==i))<10)
            temp_mask(find(temp_slice==i)) = 0;
        end
    end
    
    atlas_sumEPI_ROI(:,:,j) = temp_slice.*temp_mask;
end
% ROI time activity 
for i=1:max(atlas(:))    
    for j=1:size(fmri,4)
            tmp = fmri(:,:,:,j);
            ROI = tmp(atlas_sumEPI_ROI == i);
            sumEPI_roi_fmri(i,j) = mean(ROI(:));
    end
end
slice = 51;
figure
subplot(1,2,1),imagesc(imrotate(atlas_sumEPI(:,:,slice),90)), title('Masked Atlas with all ROIs')
subplot(1,2,2),imagesc(imrotate(atlas_sumEPI_ROI(:,:,slice),90)), title('Masked Atlas without discarded ROIs')
% Plot ROI time activity 
figure, plot(sumEPI_roi_fmri'), title('ROI time activity')
%% pt. 2 Noise regression
% load the first 6 regressors
load('./Data/FMRI/MOCOparams.mat')
% 2a regression matrix
params = [newMOCOparams WM_masked_mean_225' CSF_masked_mean_225'];
zscores = zscore(params);

%fMRI denoised
sumEPI_roi_fmri = rmmissing(sumEPI_roi_fmri,1); 
noise_regress_sumEPI =  lscov(zscores, sumEPI_roi_fmri'); 
fmri_denoised_sumEPI = sumEPI_roi_fmri' - zscores*noise_regress_sumEPI;
reg_matrix = zscores*noise_regress_sumEPI;

% plot of RC matrix before and after z-scoring
figure, 
subplot(121), imagesc(noise_regress_sumEPI), title('Regression matrix')
subplot(122), imagesc(reg_matrix), title('Regression matrix - zScore')


%% 2b  
% Cut-off 1/128
% d1 = designfilt('highpassiir','FilterOrder', 5,'PassbandFrequency',7.8125e-3);

[n,Wn] = buttord(0.0078,0.002,3,60);  
[b,a] = butter(n,Wn,'high');
filt_fmri_sumEPI =  filtfilt(b,a,fmri_denoised_sumEPI);
% lost the DC component
figure,
subplot(211), plot(fmri_denoised_sumEPI), title('ROI time activity without undesired fluctiatuions')
subplot(212), plot(filt_fmri_sumEPI), title('Signal fitered with a high-pass')
%% 3 Volyme censoring
% FD
load('./Data/FMRI/FDparams.mat')
% no volumes with a FD greater than 3.5mm
test_FD = find(FD(:,1) > 3.5);

if (length(test_FD) == 0)
    disp('no displacement greater than 3.5mm')
end
%% 4 Check pre-processing steps 
% right hippocampus region
num_ROI = 1;
figure
subplot(311), plot(sumEPI_roi_fmri(num_ROI,:)),title('Original signal for the Hippocampus'), grid on, xlim([0 225])
subplot(312), plot(fmri_denoised_sumEPI(:,num_ROI)),title('Denoised'), grid on, xlim([0 225])
% removed the slow component
subplot(313), plot(filt_fmri_sumEPI(:,num_ROI)), title('Time filtered'),grid on, xlim([0 225])

% Drift EPI
drift_fmri_sumEPI = sumEPI_roi_fmri(1,end)- sumEPI_roi_fmri(1,1);
drift_fmri_denoised_sumEPI = fmri_denoised_sumEPI(1,end)- fmri_denoised_sumEPI(1,1);
drift_fmri_time_filter_sumEPI = filt_fmri_sumEPI(1,end) - filt_fmri_sumEPI(1,1);

%% 5 FC computation
%Pairwise Pearson correlation
[conf_mat_sumEPI,P_sumEPI] =  corrcoef(filt_fmri_sumEPI);
% Functional conectivity matrix
conf_mat_sumEPI = atanh(conf_mat_sumEPI);

figure, imagesc(conf_mat_sumEPI),title('FC matrix sumEPI');
%% 6  Multiple Comparison Correction
% Using FDR
[h_valid_pVals, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P_sumEPI,0.05,'dep','yes');
%% 7 Graph Measures

% node degree
n_degree = zeros(1,size(filt_fmri_sumEPI,2));

% binarize confusion matrix 
conf_mat_bin=zeros(size(filt_fmri_sumEPI,2),size(filt_fmri_sumEPI,2));
conf_mat_bin(h_valid_pVals == 1) = 1;

color_graph = 'jet';
figure
subplot(121), imagesc(conf_mat_sumEPI), colormap(color_graph)
subplot(122), imagesc(conf_mat_bin), colormap(color_graph)
% node degree
node_deg = squeeze(sum(conf_mat_bin,1));

conf_mat_pVals = conf_mat_sumEPI.*h_valid_pVals;
% node Strength
node_strength = squeeze(sum(conf_mat_pVals,1,'omitnan'));
% node betweenness
BC=betweenness_wei(conf_mat_pVals);
%
[B_deg,I_deg] = maxk(node_deg,10);
[B_str,I_str] = maxk(node_strength,10);
[B_bet,I_bet] = maxk(BC,10);

figure
% deg
subplot(311), stem(node_deg), title('Nodes degree'), grid on, hold on, 
stem(I_deg,B_deg, 'rx'), xlim([0,60]), xticks(sort(I_deg))
% str
subplot(312), stem(node_strength), title('Node strength'), grid on, hold on,
stem(I_str,B_str, 'rx'), xlim([0,60]), xticks(sort(I_str))
% bet
subplot(313), stem(BC/max(BC)), title('Node betweenness'), grid on, hold on,
stem(I_bet,B_bet/max(BC), 'rx'),xlim([0,60]), xticks(sort(I_bet))
%% DIFFUSION MRI ANALYSIS
%% 1
% file loading 
bvals = load('./Data/DMRI/bvals');
bvecs = load('./Data/DMRI/bvecs');
diffusion_volumes_raw = load_untouch_nii('./Data/DMRI/diffusion_volumes.nii')
diffusion_brain_mask_raw = load_untouch_nii('./Data/DMRI/diffusion_brain_mask.nii')
diffusion_volumes = diffusion_volumes_raw.img;
mask = diffusion_brain_mask_raw.img;
%% 1a
nVols=size(diffusion_volumes,4);
nSlices=size(diffusion_volumes,3);
nVox=size(diffusion_volumes,1);

counter = 1;
for i=1:length(bvals)
    if(bvals(i) ~=0)
        new_bvals_no_0(counter) = bvals(i);
        counter = counter + 1;
    end
end

shells_counter = 0;

shells(1) = new_bvals_no_0(1);
shells_counter = shells_counter + 1;

for i=2:length(new_bvals_no_0)
    tmp = new_bvals_no_0(i);
    guard = 1;
    for k=1:shells_counter
        if (tmp < (shells(k) + 20) && tmp > (shells(k) - 20))
            guard = 0;
        end
    end
    
    if(guard == 1)
        shells(i) = tmp;
        shells_counter = shells_counter + 1;
    end
end
DWI_number = size(bvals,2);
% results: 
% DWI = 103
% shells_counter = 2

%% 1.b (*)
% selected voxel for which the mean FA value is close to 0: (81,82) in
% slice 50
% FA close to zero = isotropic movement of water molecules, so CSF

diff_signal = squeeze(diffusion_volumes(81,82,50,:));

%% 1.c inter/intra b-value variabilities
figure, 
subplot(211), plot(diff_signal), title('Mean diffusion signal for a single voxel'), xlim([0 103])
% sort b-values
[B_sort_tmp, I_sort_tmp] = sort(bvals);
diff_signal = diff_signal(I_sort_tmp);
subplot(212), plot(diff_signal), title('Mean diffusion signal for a single voxel - sorted'), xlim([0 103])
%% 2a 
% find closest bval to 1000
min_dist_1000_bval = 0;
min_dist_1000 = 1000; 
for i=1:length(bvals)
    dist = abs(bvals(i) - 1000); 
    if(dist < min_dist_1000)
        min_dist_1000 = dist;
        min_dist_1000_bval = bvals(i);
    end
end
% it's about 700

% Create the new matrix
counter = 1;
counter_minDist = 1;
counter_onlyZero = 1;
for i=1:nVols
    bvals_tmp = bvals(i);
    if(bvals_tmp == 0 | (bvals_tmp < (min_dist_1000_bval + 20) && bvals_tmp > (min_dist_1000_bval  - 20)))
        new_diffusion_matrix(:,:,:,counter) = diffusion_volumes(:,:,:,i);
        new_bvals(counter) = bvals(i);
        new_bvecs(:,counter) = bvecs(:,i);
        counter = counter + 1;
    end
    % Matrix with only bvals = 0
    if(bvals_tmp == 0)
        new_diffusion_matrix_onlyZero(:,:,:,counter_onlyZero) = diffusion_volumes(:,:,:,i);
        new_bvals_onlyZero(counter_onlyZero) = bvals(i);
        new_bvecs_onlyZero(:,counter_onlyZero) = bvecs(:,i);
        counter_onlyZero = counter_onlyZero + 1;
    end
end
%% 2b

% Build the B design matrix for the linear least squares approach
for ii=1:length(new_bvals)
    
    B(ii,1)=new_bvecs(1,ii)^2;
    B(ii,2)=new_bvecs(2,ii)^2;
    B(ii,3)=new_bvecs(3,ii)^2;
    B(ii,4)=new_bvecs(1,ii)*new_bvecs(2,ii);
    B(ii,5)=new_bvecs(1,ii)*new_bvecs(3,ii);
    B(ii,6)=new_bvecs(2,ii)*new_bvecs(3,ii);
end


B=new_bvals'.*B;

% S0 is the voxel wise value of the first b=0 volume
S0 = new_diffusion_matrix_onlyZero(:,:,:,1);

% Normalize the signals and pass to the logarithm for the fit
for ii=1:length(new_bvals)
    Slog(:,:,:,ii)=log((new_diffusion_matrix(:,:,:,ii)./(S0))+eps);
end

nVols=size(new_diffusion_matrix,4);
nSlices=size(new_diffusion_matrix,3);
nVox=size(new_diffusion_matrix,1);


%initialize the structures which will be used to contain DTI parameters
FA=zeros(nVox,nVox,nSlices);
MD=zeros(nVox,nVox,nSlices);

FirstX=zeros(nVox,nVox,nSlices);
FirstY=zeros(nVox,nVox,nSlices);
FirstZ=zeros(nVox,nVox,nSlices);
%%
% start the cycle to fit the voxel-wise diffusion tensor
for jj=1:nSlices
    
    %print fitting progress
    disp([' Fitting Slice ',num2str(jj)])
    
    for kk=1:1:nVox
        for ll=1:1:nVox
            
            %check if current voxel belongs to the mask or has 0 info to
            %give

            
            if (mask(kk,ll,jj) && S0(kk,ll,jj)~=0)
                
                %extract the signal from each voxel
                VoxSignal=squeeze(Slog(kk,ll,jj,:));
                
                %fit the DTI
                D=-inv(B'*B)*B'*VoxSignal;
                
                %reconstruct the diffusion tensor from the fitted
                %parameters
                T=[D(1) D(4)/2 D(5)/2;
                    D(4)/2 D(2) D(6)/2;
                    D(5)/2 D(6)/2 D(3)];
                
                %T(isnan(T))=0;
                
                %compute eigenvalues and eigenvectors
                [eigenvects, eigenvals]=eig(T);
                eigenvals=diag(eigenvals);
                
                
                %Manage negative eigenvals
                % if all <0 -> take the absolute value
                %otherwise -> put negatives to zero
                if((eigenvals(1)<0)&&(eigenvals(2)<0)&&(eigenvals(3)<0)), eigenvals=abs(eigenvals);end
                if(eigenvals(1)<0), eigenvals(1)=0; end
                if(eigenvals(2)<0), eigenvals(2)=0; end
                if(eigenvals(3)<0), eigenvals(3)=0; end
                
                First = eigenvects(:,3);
                
                %compute FA
                FAv=(1/sqrt(2))*( sqrt((eigenvals(1)-eigenvals(2)).^2+(eigenvals(2)-eigenvals(3)).^2 + ...
                    (eigenvals(1)-eigenvals(3)).^2)./sqrt(eigenvals(1).^2+eigenvals(2).^2+eigenvals(3).^2) );
                FA(kk,ll,jj)=FAv;
                
                %Compute the MD
                MDv=(eigenvals(1)+eigenvals(2)+eigenvals(3))/3;
                MD(kk,ll,jj)=MDv;
                
                %sort eigenvalues, eigenvectors
                [sorted,idx_sort]=sort(eigenvals);
                
                eigenvals=eigenvals(idx_sort);
                eigenvects=eigenvects(:,idx_sort);
                
                %take principal eigenvector, decompose it
                First = eigenvects(:,3);
                
                FirstX(kk,ll,jj)=abs(First(1))*FAv;
                FirstY(kk,ll,jj)=abs(First(2))*FAv;
                FirstZ(kk,ll,jj)=abs(First(3))*FAv;
                   
            end
        end
    end
end
FA_1 = FA;
MD_1 = MD;

%%
clear color
% create the colour-encoded directional map
color(:,:,:,1)=FirstX;
color(:,:,:,2)=FirstY;
color(:,:,:,3)=FirstZ;
figure
%visualize the map
for ii=1:1:50
    slice=reshape(color(:,:,ii,:),nVox,nVox,3);
    image(imrotate(slice,90)),title('Map visualization for 2.b')
    pause(0.1)
    
end
%% 2.b
clear Slog
clear S0
% S0 is the voxel wise mean value of all the b=0 volumes
S0 = mean(new_diffusion_matrix_onlyZero,4);

%normalize the signals and pass to the logarithm for the fit
for ii=1:length(new_bvals)
    Slog(:,:,:,ii)=log((new_diffusion_matrix(:,:,:,ii)./(S0))+eps);
end

nVols=size(new_diffusion_matrix,4);
nSlices=size(new_diffusion_matrix,3);
nVox=size(new_diffusion_matrix,1);


%initialize the structures which will be used to contain DTI parameters
FA=zeros(nVox,nVox,nSlices);
MD=zeros(nVox,nVox,nSlices);

FirstX=zeros(nVox,nVox,nSlices);
FirstY=zeros(nVox,nVox,nSlices);
FirstZ=zeros(nVox,nVox,nSlices);
% start the cycle to fit the voxel-wise diffusion tensor
for jj=1:nSlices
    
    %print fitting progress
    disp([' Fitting Slice ',num2str(jj)])
    
    for kk=1:1:nVox
        for ll=1:1:nVox
            
            %check if current voxel belongs to the mask or has 0 info to
            %give
            if (mask(kk,ll,jj) && S0(kk,ll,jj)~=0)
                
                %extract the signal from each voxel
                VoxSignal=squeeze(Slog(kk,ll,jj,:));
                
                %fit the DTI
                D=-inv(B'*B)*B'*VoxSignal;
                
                %reconstruct the diffusion tensor from the fitted
                %parameters
                T=[D(1) D(4)/2 D(5)/2;
                    D(4)/2 D(2) D(6)/2;
                    D(5)/2 D(6)/2 D(3)];
                
                %T(isnan(T))=0;
                
                %compute eigenvalues and eigenvectors
                [eigenvects, eigenvals]=eig(T);
                eigenvals=diag(eigenvals);
                
                
                %Manage negative eigenvals
                % if all <0 -> take the absolute value
                %otherwise -> put negatives to zero
                if((eigenvals(1)<0)&&(eigenvals(2)<0)&&(eigenvals(3)<0)), eigenvals=abs(eigenvals);end
                if(eigenvals(1)<0), eigenvals(1)=0; end
                if(eigenvals(2)<0), eigenvals(2)=0; end
                if(eigenvals(3)<0), eigenvals(3)=0; end
                
                First = eigenvects(:,3);
                
                %compute FA
                FAv=(1/sqrt(2))*( sqrt((eigenvals(1)-eigenvals(2)).^2+(eigenvals(2)-eigenvals(3)).^2 + ...
                    (eigenvals(1)-eigenvals(3)).^2)./sqrt(eigenvals(1).^2+eigenvals(2).^2+eigenvals(3).^2) );
                FA(kk,ll,jj)=FAv;
                
                %Compute the MD
                MDv=(eigenvals(1)+eigenvals(2)+eigenvals(3))/3;
                MD(kk,ll,jj)=MDv;
                
                %sort eigenvalues, eigenvectors
                [sorted,idx_sort]=sort(eigenvals);
                
                eigenvals=eigenvals(idx_sort);
                eigenvects=eigenvects(:,idx_sort);
                
                %take principal eigenvector, decompose it
                First = eigenvects(:,3);
                
                FirstX(kk,ll,jj)=abs(First(1))*FAv;
                FirstY(kk,ll,jj)=abs(First(2))*FAv;
                FirstZ(kk,ll,jj)=abs(First(3))*FAv;
                   
            end
        end
    end
end

FA_2 = FA;
MD_2 = MD;

%%
clear color
% create the colour-encoded directional map
color(:,:,:,1)=FirstX;
color(:,:,:,2)=FirstY;
color(:,:,:,3)=FirstZ;

figure
%visualize the map
for ii=1:1:50
    slice=reshape(color(:,:,ii,:),nVox,nVox,3);
    image(imrotate(slice,90)),title('Map visualization for 2.c')
    pause(0.1)
end

%% 2.d Voxel-wise coefficients of variation

mean_FA_1 = mean2(FA_1);
mean_MD_1 = mean2(MD_1);
mean_FA_2 = mean2(FA_2);
mean_MD_2 = mean2(MD_2);

% VC for FAs and MDs
for i=1:size(FA_1,3)
    for j=1:size(FA_1,1)
        for k=1:size(FA_1,2)
            CV_FA_1(j,k,i) = 100 * (mean_FA_1 - FA_1(j,k,i)) ./ mean_FA_1;
            CV_MD_1(j,k,i) = 100 * (mean_MD_1 - MD_1(j,k,i)) ./ mean_MD_1;
            CV_FA_2(j,k,i) = 100 * (mean_FA_2 - FA_2(j,k,i)) ./ mean_FA_2;
            CV_MD_2(j,k,i) = 100 * (mean_MD_2 - MD_2(j,k,i)) ./ mean_MD_2;
        end
    end
end

%for i=1:size(CV_FA_1,3)
slice = 72;
colorbar_color = 'jet';
figure
subplot(221), imagesc(imrotate(CV_FA_1(:,:,slice),90)),
colormap(colorbar_color),colorbar, title('CV FA 1')
subplot(222), imagesc(imrotate(CV_FA_2(:,:,slice),90)),
colormap(colorbar_color),colorbar, title('CV FA 2')
subplot(223), imagesc(imrotate(CV_MD_1(:,:,slice),90)),
colormap(colorbar_color),colorbar, title('CV MD 1')
subplot(224), imagesc(imrotate(CV_MD_2(:,:,slice),90)),
colormap(colorbar_color),colorbar, title('CV MD 2')
%pause
%end
%
colorbar_color = 'jet';
figure
subplot(221), imagesc(imrotate(abs(mean(CV_FA_1-CV_FA_2,3)),90)),
colormap(colorbar_color),colorbar, title('CV FA abs difference on mean')
subplot(223), imagesc(imrotate(abs(mean((CV_MD_1-CV_MD_2),3)),90)),
colormap(colorbar_color),colorbar, title('CV MD abs difference on mean')

subplot(222), imagesc(imrotate(abs(CV_FA_1(:,:,slice)-CV_FA_2(:,:,slice)),90)),
colormap(colorbar_color),colorbar, title('CV FA abs difference single slice')
subplot(224), imagesc(imrotate(abs(CV_MD_1(:,:,slice)-CV_MD_2(:,:,slice)),90)),
colormap(colorbar_color),colorbar, title('CV MD abs difference single slice')
% Looking at FA and MD difference, FA is the most affectedby the
% different normalization choise. This may vary looking at
% different slices 
%%  2.e
% FA_2, MD_2 as S0 estimated as a mean should be better. 

% Plot the selected FA and MD for a slice
figure
subplot(121), imagesc(imrotate(FA_2(:,:,slice),90),[0, 1]),colormap(colorbar_color),colorbar, title('FA slice 72')
subplot(122), imagesc(imrotate(MD_2(:,:,slice),90),[0, 5e-3]),colormap(colorbar_color),colorbar, title('MD slice 72')

% Saving maps
save_3D_nii(FA_2,'./Data/DMRI/diffusion_brain_mask.nii' ,'./Saved_data/FA.nii');
save_3D_nii(MD_2,'./Data/DMRI/diffusion_brain_mask.nii' ,'./Saved_data/MD.nii');
%% 2.f
tmp = 50;
% mask FA and MD with GM and sumEPI
masked_FA = FA_2.*GM_erode;
masked_FA = masked_FA.*sumEPI_mask;
masked_MD = MD_2.*GM_erode;
masked_MD = masked_MD.*sumEPI_mask;

% visualize reults
figure, 
subplot(221),imagesc(imrotate(FA_2(:,:,tmp),90),[0, 1]),colormap(colorbar_color),colorbar,title('FA')
subplot(223),imagesc(imrotate(masked_FA(:,:,tmp),90),[0 1]),colormap(colorbar_color),colorbar, title('masked FA')
subplot(222),imagesc(imrotate(MD_2(:,:,tmp),90),[0 5e-3]),colormap(colorbar_color),colorbar,title('MD')
subplot(224),imagesc(imrotate(masked_MD(:,:,tmp),90),[0 5e-3]),colormap(colorbar_color),colorbar,title('masked MD')

% Mean values extraction in each ROI for FA and MD maps 
for i=1:max(atlas_sumEPI(:))    
        %FA    
        ROI = masked_FA(atlas_sumEPI_ROI == i);
        sumEPI_roi_FA(i) = mean(ROI(:));
        %MD
        ROI = masked_MD(atlas_sumEPI_ROI == i);
        sumEPI_roi_MD(i) = mean(ROI(:));
end


%% back toi point 1b (*)
% To select a voxel in FA populated mainly with CSF

figure, imagesc(FA_2(:,:,50));
% Selected (81,82) of slice 50


%% DMRI/fMRI integration
x_ax = 1:1:length(node_deg);
% remove NaN from FA and MD
FA_ROI_2 = rmmissing(sumEPI_roi_FA,2);
MD_ROI_2 = rmmissing(sumEPI_roi_MD,2);
node_between_norm = (BC/max(BC))';

% Visual inspection of scatter plots
figure
subplot(231), scatter(FA_ROI_2, node_deg), hold on, title('Node degree and FA'), grid on
subplot(232), scatter(FA_ROI_2, node_strength), title('Node strength and FA'), grid on
subplot(233), scatter(FA_ROI_2, node_between_norm), title('Node betweenness and FA'), grid on
subplot(234), scatter(MD_ROI_2, node_deg), title('Node degree and MD'), grid on
subplot(235), scatter(MD_ROI_2, node_strength), title('Node strength and MD'), grid on
subplot(236), scatter(MD_ROI_2, node_between_norm), title('Node betweenness and MD'), grid on

% Pearson Correlation
rho_deg_FA = corrcoef(node_deg,FA_ROI_2);
rho_streng_FA = corrcoef(node_strength,FA_ROI_2);
rho_between_FA = corrcoef(node_between_norm,FA_ROI_2);
rho_deg_MD = corrcoef(node_deg,MD_ROI_2);
rho_streng_MD = corrcoef(node_strength,MD_ROI_2);
rho_between_MD = corrcoef(node_between_norm,MD_ROI_2);  
%% Save the required files 
GC = reg_matrix;
filtered_regression_and_temporal_fmri = filt_fmri_sumEPI;
zFC = conf_mat_sumEPI;
zFC_corr = conf_mat_pVals;
DEG = node_deg;
STR = node_strength;
BTW_NORM = node_between_norm;
cv_MD = CV_MD_2;
cv_FA = CV_FA_2;
save('gitMelk_variables.mat', 'GC','filtered_regression_and_temporal_fmri', ...
        'zFC','zFC_corr','DEG','STR','BTW_NORM','cv_MD','cv_FA')
