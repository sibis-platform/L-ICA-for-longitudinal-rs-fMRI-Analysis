function singleSubjectMultiICA_general(subjectIdx,filenames, mask, IC, age, outputPath, Options)
%singleSubjectMultiICA - perform subject-specific longitudinal ICA 
%
% Inputs:
%    filenames - A cell array of BOLD signal (pre-processed NIFTI) filenames.
%    Currently, we only support N (array length) == 3, i.e., 3 BOLD images
%    for 3 subject-specific longitudinal visits.
%
%    mask - binary NIFTI
%
%    IC - group-level ICs (4D NIFTI)
%
%    age - age of the subject associated with the sessions
%  
%    outputPath - path for output subject-specific IFN files
%
% Optional Inputs:
%
%    ICAOption - 1 (RL-ICA), 2 (VL-ICA, default), 3 (Both)
%
%    parcellation - for RL-ICA, a parcellation of the brain must be
%    provided; the parcel indices are assumed to be 1:N for N parcels
%
%    icIdx - an 1D array for indicating which of the group-level ICs correspond
%    to real IFN.
%
%    k1, k2 - weighting parameters for RL-ICA (k1) or VL-ICA (k1,k2)
%
% To do:
%    numEigs: dimension reduction for individual BOLD, the reduced dimension
%    should adapt to specific datasets, e.g., retaining 95% energy, 99%
%    energy, etc.
%
%    
% Example: 
%    filenames = {'a.nii.gz','b.nii.gz','c.nii.gz'};
%    age = [1,2,3];
%    mask = 'mask.nii.gz';
%    IC = 'melodic.nii.gz';  % FSL melodic output, assume it contains 5 group-level ICs
%    outputPath = './results';
%    Options.ICAOption = 3;
%    Options.parcellation = '100ROI.nii.gz';
%    Options.icIdx = [1,2,5];   % 3 of the 5 ICs correspond to true IFNs
%    singleSubjectMultiICA(filenames, mask, IC, age, outputPath, Options);
%
% Other m-files required: ./tools/*.m
%
% Author: Qingyu Zhao, Ph.D., 
% Stanford University / SRI international
% email address: qingyuz@stanford.edu
% Last Modified: May, 2018
    
    %% load/initialize input parameters
    icIdx = 1:size(IC,4);
    ICAOption = 2;
    k1 = 0.2;
    k2 = 100;
    
    if isfield(Options, 'icIdx')
        icIdx = Options.icIdx;
    end
	if isfield(Options, 'ICAOption')
        ICAOption = Options.ICAOption;
    end
    if isfield(Options, 'k1')
        k1 = Options.k1;
    end
	if isfield(Options, 'k2')
        k2 = Options.k2;
    end
    if isfield(Options, 'parcellation')
        parcellation = load_nii(Options.parcellation);
        parcelNum = max(parcellation.img(:));
    end
        
    maskVoxNum = sum(mask.img(:) > 0);
    ICNum = size(IC.img,4);

    regression_x = zeros(ICNum,maskVoxNum);
    for ic = 1:ICNum
        ic_img = squeeze(IC.img(:,:,:,ic));
        regression_x(ic,:) = ic_img(mask.img > 0)';
    end
    
    %% load BOLD signals
    for i = 1:length(filenames)
        raw{i} = load_nii(filenames{i});
    end
    sessionNum = length(filenames);
    
    %% generate RL-ICA constraint MASK   
    if (ICAOption == 1) || (ICAOption == 3)
        S = zeros(parcelNum,sum(mask.img(:) > 0));
        p = parcellation.img(mask.img > 0);

        for i = 1:parcelNum
            idx = (p == i);

            S(i,:) = idx;
            S(i,:) = S(i,:) / sum(S(i,:));
        end 
    end
    
    %% ICA preparation: subject-specific PCA, rescale group-level templates to adapt to subject-specific data
    numEigs = 80; 
    
    % numEigs: dimension reduction for individual BOLD, the reduced dimension
    % should adapt to specific datasets, e.g., retaining 95% energy, 99%
    % energy, etc.
    for i = 1:sessionNum
        timeCourseLen{i} = size(raw{i}.img,4);
        rawData = zeros(timeCourseLen{i}, maskVoxNum);
        for t = 1:timeCourseLen{i}
            raw_img = squeeze(raw{i}.img(:,:,:,t));
            rawData(t,:) = raw_img(mask.img > 0)';
        end    
        [whitenedData{i},w_init{i},ref{i},numSamples] = whiteningDataMulti(mask,rawData,regression_x,icIdx,numEigs);
    end
    
    % to use a single group guidance for all 
    refAvg = ref{1};
    for i = 2:sessionNum
        refAvg = refAvg + ref{i};
    end
    refAvg = refAvg / sessionNum;
    
    %% RL-ICA
    if (ICAOption == 1) || (ICAOption == 3)
        w = RLICA(whitenedData, refAvg, S, w_init, age, k1);

        for i = 1:length(icIdx)
            s = zeros(maskVoxNum,1);
            for j = 1:sessionNum
                si{j} = whitenedData{j}'*w{j}(:,i);
                s = s + si{j};
            end
            s = s/sessionNum;
            m = mean(s);
            v = std(s);
            
            s = [];
            for j = 1:sessionNum
                %si{j} = (si{j}-m)/v;
                si{j} = si{j}/norm(si{j});
                s = [s,si{j}];
                filename = sprintf('%s/IC_%d_%d_%s_R.nii.gz',outputPath,icIdx(i),j,subjectIdx);
                new_ic = mask;
                new_ic.img(mask.img>0) = si{j};
                save_nii(new_ic,filename);
            end
    
            b = zeros(1,maskVoxNum);
            for bidx = 1:length(s)
                beta = [ones(sessionNum,1),age]\ s(bidx,:)';
                b(bidx) = beta(2);
            end
            
            slope = mask;
            slope.img(mask.img>0) = b;
            filename = sprintf('%s/slope_%d_%s.nii_R.nii.gz',outputPath,icIdx(i),subjectIdx);
            save_nii(slope,filename);
        end
    end
    
    %% VL-ICA
    if (ICAOption == 2) || (ICAOption == 3)
        w = VLICA(whitenedData, refAvg, w_init, age, k1, k2);

        for i = 1:length(icIdx)
            s = zeros(maskVoxNum,1);
            for j = 1:sessionNum
                si{j} = whitenedData{j}'*w{j}(:,i);
                s = s + si{j};
            end
            s = s/sessionNum;
            m = mean(s);
            v = std(s);
            
            s = [];
            for j = 1:sessionNum
                %si{j} = (si{j}-m)/v;
                si{j} = si{j}/norm(si{j});
                s = [s,si{j}];
                filename = sprintf('%s/IC_%d_%d_%s_V_%f_%f.nii.gz',outputPath,icIdx(i),j,subjectIdx,k1,k2);
                new_ic = mask;
                new_ic.img(mask.img>0) = si{j};
                save_nii(new_ic,filename);
            end
            
            b = zeros(1,maskVoxNum);
            for bidx = 1:length(s)
                beta = [ones(sessionNum,1),age]\ s(bidx,:)';
                b(bidx) = beta(2);
            end
            
            slope = mask;
            slope.img(mask.img>0) = b;
            filename = sprintf('%s/slope_%d_%s.nii_V_%f_%f.nii.gz',outputPath,icIdx(i),subjectIdx,k1,k2);
            save_nii(slope,filename);
        end

    end
end
