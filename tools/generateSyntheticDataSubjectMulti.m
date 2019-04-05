clear;

loadedData = false;
regenerateTemplate = true;
ICAOption = 3;
ICARegularization = [0.05,1];    
icIdx = [2,3,5,6,11,18];

if ~loadedData
addpath('FastICA_25');
addpath('/fs/neurosci01/qingyuz/rsfmri/Seed2Vox/nifti');

mask = load_nii('/fs/neurosci01/djk/ncanda/atlas/sri24/restingstate/baseline/melodic/results_6mm_d25_dx/mask.nii.gz');
IC = load_nii('/fs/neurosci01/djk/ncanda/atlas/sri24/restingstate/baseline/melodic/results_6mm_d25_dx/melodic_IC.nii.gz');
parcellation = load_nii('masks/sri24_functional_parcellation_100.nii.gz');

maskVoxNum = sum(mask.img(:) > 0);
regression_x = zeros(25,maskVoxNum);
for ic = 1:25
	ic_img = squeeze(IC.img(:,:,:,ic));
	regression_x(ic,:) = ic_img(mask.img > 0)';
end
end

data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline_ex0.csv');
data_baseline_selected = data_baseline((strcmp(data_baseline.dx,'ctrl') | ...
                                        strcmp(data_baseline.dx,'heavy'))&...
                                       data_baseline.b_restingstate == 1, :);
                                   
data_f1y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_1y_ex0.csv');
data_f1y_selected = data_f1y((strcmp(data_f1y.dx,'ctrl') | ...
                                          strcmp(data_f1y.dx,'heavy'))&...
                                          data_f1y.b_restingstate == 1, :);
                                      
                                      
data_f2y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_2y_ex0.csv');
data_f2y_selected = data_f2y((strcmp(data_f2y.dx,'ctrl') | ...
                                          strcmp(data_f2y.dx,'heavy'))&...
                                          data_f2y.b_restingstate == 1, :);

subjects_02 = intersect(data_baseline_selected.subject,data_f2y_selected.subject);
subjects_01 = intersect(data_baseline_selected.subject,data_f1y_selected.subject);
subjects_012 = intersect(data_f2y_selected.subject,subjects_01);

%% generate constraint MASK
XC = zeros(size(mask.img,1),size(mask.img,1),size(mask.img,1));
YC = zeros(size(mask.img,1),size(mask.img,1),size(mask.img,1));
ZC = zeros(size(mask.img,1),size(mask.img,1),size(mask.img,1));

for i = 1:size(mask.img,1)
    for j = 1:size(mask.img,2)
        for k = 1:size(mask.img,3)
            XC(i,j,k) = i;
            YC(i,j,k) = j;
            ZC(i,j,k) = k;
        end
    end
end

parcelNum = max(parcellation.img(:));
S = zeros(parcelNum,sum(mask.img(:) > 0));
p = parcellation.img(mask.img > 0);
xc = XC(mask.img > 0);
yc = YC(mask.img > 0);
zc = ZC(mask.img > 0);

for i = 1:max(parcellation.img(:))
    idx = (p == i);
    S(i,idx) = 1;
    mxc = mean(xc(idx));
    myc = mean(yc(idx));
    mzc = mean(zc(idx));
    
    dis = (xc - mxc).^2 + (yc - myc).^2 + (zc - mzc).^2;
    region = (dis <= 300);
    
    S(i,:) = exp(-dis / 50);
    S(i,:) = S(i,:) / norm(S(i,:));

end

%% Begin Experiments
for subjectIdx = 1:length(subjects_012)
    %% randomly select ROI region to modify
    if regenerateTemplate
        roiSelected = zeros(100,2);
        for i = 1:25
            randROI = randperm(100);
            roiSelected(i,1) = randROI(1);
            roiSelected(i,2) = randROI(2);
        end
        filename = sprintf('roiSelected_%s.mat',subjects_012{subjectIdx});
        save(filename,'roiSelected');
    end
    
    filename1 = sprintf('/fs/ncanda-share/cases/%s/standard/baseline/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
	filename2 = sprintf('/fs/ncanda-share/cases/%s/standard/followup_1y/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
	filename3 = sprintf('/fs/ncanda-share/cases/%s/standard/followup_2y/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
    if ~loadedData
        raw1 = load_nii(filename1);
        raw2 = load_nii(filename2);    
        raw3 = load_nii(filename3);
    end
    
%%  create subject template
	if regenerateTemplate
    rawData = zeros(size(raw1.img,4)+size(raw2.img,4)+size(raw3.img,4),maskVoxNum);
    for t = 1:size(raw1.img,4)
        raw_img = squeeze(raw1.img(:,:,:,t));
        rawData(t,:) = raw_img(mask.img > 0)';
    end
    for t = 1:size(raw2.img,4)
        raw_img = squeeze(raw2.img(:,:,:,t));
        rawData(t+size(raw1.img,4),:) = raw_img(mask.img > 0)';
    end
    for t = 1:size(raw3.img,4)
        raw_img = squeeze(raw3.img(:,:,:,t));
        rawData(t+size(raw1.img,4)+size(raw2.img,4),:) = raw_img(mask.img > 0)';
    end   
    w_init = regression_x'\rawData';
    regression_x_subject = w_init'\rawData;
    rawData1 = rawData(1:size(raw1.img,4),:);
    rawData2 = rawData(1+size(raw1.img,4):size(raw1.img,4)+size(raw2.img,4),:);
    rawData3 = rawData(1+size(raw1.img,4)+size(raw2.img,4):end,:);
    
%     parcelIdx = parcellation.img(mask.img > 0);
%     regression_x_subject_2 = regression_x_subject;
%     for i = 1:25
%         regression_x_subject_2(i,parcelIdx == roiSelected(i,1)) = regression_x_subject_2(i, parcelIdx == roiSelected(i,1)) + 0.3;
%         regression_x_subject_2(i,parcelIdx == roiSelected(i,2)) = regression_x_subject_2(i, parcelIdx == roiSelected(i,2)) - 0.3;
%     end
%     
%     regression_x_subject_3 = regression_x_subject;
%     for i = 1:25
%         regression_x_subject_3(i,parcelIdx == roiSelected(i,1)) = regression_x_subject_3(i, parcelIdx == roiSelected(i,1)) + 0.6;
%         regression_x_subject_3(i,parcelIdx == roiSelected(i,2)) = regression_x_subject_3(i, parcelIdx == roiSelected(i,2)) - 0.6;
%     end
%     
%     w_init1 = w_init(:,1:size(raw1.img,4));
%     w_init2 = w_init(:,1+size(raw1.img,4):size(raw1.img,4)+size(raw2.img,4));
%     w_init3 = w_init(:,1+size(raw1.img,4)+size(raw2.img,4):end);
% 
%     rawData1 = w_init1' * regression_x_subject   + randn(size(raw2.img,4),maskVoxNum)*4;%rawData1 - w_init1' * regression_x_subject;% 
%     rawData2 = w_init2' * regression_x_subject_2 + randn(size(raw2.img,4),maskVoxNum)*4;%rawData2 - w_init2' * regression_x_subject;% randn(size(raw2.img,4),maskVoxNum)*4;
%     rawData3 = w_init3' * regression_x_subject_3 + randn(size(raw2.img,4),maskVoxNum)*4;%rawData3 - w_init3' * regression_x_subject;% randn(size(raw3.img,4),maskVoxNum)*4;
%     rawData = [rawData1;rawData2;rawData3];
%     
%     for i = 1:length(icIdx)
%         filename = sprintf('IC_%d_1g_%s.nii.gz',icIdx(i),subjects_012{subjectIdx});
%         new_ic = mask;
%         new_ic.img(mask.img>0) = regression_x_subject(icIdx(i),:);
%         save_nii(new_ic,filename);
%         filename = sprintf('IC_%d_2g_%s.nii.gz',icIdx(i),subjects_012{subjectIdx});
%         new_ic = mask;
%         new_ic.img(mask.img>0) = regression_x_subject_2(icIdx(i),:);
%         save_nii(new_ic,filename);
%         filename = sprintf('IC_%d_3g_%s.nii.gz',icIdx(i),subjects_012{subjectIdx});
%         new_ic = mask;
%         new_ic.img(mask.img>0) = regression_x_subject_3(icIdx(i),:);
%         save_nii(new_ic,filename);
%     end
	end
    
%%  ICA preparation
    numEigs = 80;
    
    [whitenedData1,w_init1,ref1,numSamples] = whiteningDataMulti(mask,rawData1,regression_x,icIdx,numEigs);
    [whitenedData2,w_init2,ref2,numSamples] = whiteningDataMulti(mask,rawData2,regression_x,icIdx,numEigs);
    [whitenedData3,w_init3,ref3,numSamples] = whiteningDataMulti(mask,rawData3,regression_x,icIdx,numEigs);
    
    ref = (ref1 + ref2 + ref3)/3;
    ref1 = ref; ref2 = ref; ref3 = ref;
    
    X1 = whitenedData1; X2 = whitenedData2; X3 = whitenedData3;
    XX1 = X1*X1'; XX2 = X2*X2'; XX3 = X3*X3';
    X1X2 = X1*X2'; X1X3 = X1*X3'; X2X3 = X2*X3';
	A = S*[-X1',2*X2',-X3'];
    v = randn(numSamples,1); v = (v - mean(v))/std(v);
    v = repmat(v,1,length(icIdx));
	exv = exp(-v.^2/2);
%% ICA longitudinal constraint
    if (ICAOption == 1) || (ICAOption == 3)

    ro2 = 1 / (numEigs *3);    
    for k = ICARegularization
        ro1 = max(4 / parcelNum, 1 / parcelNum * k);
        alpha = k/numSamples;

        w1 = w_init1; 
        w2 = w_init2; 
        w3 = w_init3; 
        
        p1 = w1; p2 = w2; p3 = w3; b1 = zeros(parcelNum,length(icIdx)); b2 = zeros(numEigs*3,length(icIdx));
        iter = 1;
        while iter < 12
                options = struct('GradObj','on','Display','final','LargeScale','on','HessUpdate','lbfgs','InitialHessType','identity','MaxIter',200);
                %fun = @(w)primalProblem(w,X1,X2,X3,XX1,XX2,XX3,X1X2,X2X3,X1X3,ref1,ref2,ref3,v,alpha,nu,lambda1,lambda2,lambda3,lambda,numSamples);
                p = [p1;p2;p3];
                AA = A'*A;
                Ab1 = A'*b1;

                fun = @(w)primalProblemUncMulti(w,X1,X2,X3,XX1,XX2,XX3,A,b1,AA,Ab1,p,b2,ref1,ref2,ref3,exv,alpha,ro1,ro2,numSamples,numEigs);
                w = fminlbfgs(fun,[p1;p2;p3],options);
                w1 = w(1:numEigs,:); w2 = w(numEigs+1:numEigs*2,:); w3 = w(numEigs*2+1:numEigs*3,:);

                [u,s,v] = svd(w1+b2(1:numEigs,:));                p1 = u * eye(numEigs,length(icIdx)) * v'; 
                [u,s,v] = svd(w2+b2(numEigs+1:numEigs*2,:));      p2 = u * eye(numEigs,length(icIdx)) * v'; 
                [u,s,v] = svd(w3+b2(numEigs*2+1:numEigs*3,:));    p3 = u * eye(numEigs,length(icIdx)) * v'; 
                
                b1 = b1 + A*w;
                b2 = b2 + (w - [p1;p2;p3]);
                ro1 = min(ro1 * 1,60/parcelNum * max(k/2,1));
                ro2 = min(ro2 * 1,20/(numEigs *3) * max(k/2,1));
                
                c1 = w1'*w1 - eye(length(icIdx)); 
                c2 = w2'*w2 - eye(length(icIdx)); 
                c3 = w3'*w3 - eye(length(icIdx));
                c4 = sum(sum(abs(A*w)))/parcelNum;
                fprintf('%d:%f: c1:%f c2:%f c3:%f c4:%f\n',iter,k, sum(sum(abs(c1))),sum(sum(abs(c2))),sum(sum(abs(c3))),c4);

                iter = iter + 1;
        end
    
        for i = 1:length(icIdx)
            filename = sprintf('./results_multi/IC_%d_1_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X1'*w1(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename);
            filename = sprintf('./results_multi/IC_%d_2_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X2'*w2(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename);
            filename = sprintf('./results_multi/IC_%d_3_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X3'*w3(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename); 
        end
    end
    end
    
%% ICA Independent  
    if (ICAOption == 2) || (ICAOption == 3)
    ro2 = 1 / (numEigs *3);    
    for k = ICARegularization
        ro1 = 0;
        alpha = k/numSamples;

        w1 = w_init1; 
        w2 = w_init2; 
        w3 = w_init3; 
        
        p1 = w1; p2 = w2; p3 = w3; b1 = zeros(parcelNum,length(icIdx)); b2 = zeros(numEigs*3,length(icIdx));
        iter = 1;
        while iter < 10
                options = struct('GradObj','on','Display','final','LargeScale','on','HessUpdate','lbfgs','InitialHessType','identity','MaxIter',200);
                %fun = @(w)primalProblem(w,X1,X2,X3,XX1,XX2,XX3,X1X2,X2X3,X1X3,ref1,ref2,ref3,v,alpha,nu,lambda1,lambda2,lambda3,lambda,numSamples);
                p = [p1;p2;p3];
                AA = A'*A;
                Ab1 = A'*b1;

                fun = @(w)primalProblemUncMulti(w,X1,X2,X3,XX1,XX2,XX3,A,b1,AA,Ab1,p,b2,ref1,ref2,ref3,exv,alpha,ro1,ro2,numSamples,numEigs);
                w = fminlbfgs(fun,[p1;p2;p3],options);
                w1 = w(1:numEigs,:); w2 = w(numEigs+1:numEigs*2,:); w3 = w(numEigs*2+1:numEigs*3,:);

                [u,s,v] = svd(w1+b2(1:numEigs,:));                p1 = u * eye(numEigs,length(icIdx)) * v'; 
                [u,s,v] = svd(w2+b2(numEigs+1:numEigs*2,:));      p2 = u * eye(numEigs,length(icIdx)) * v'; 
                [u,s,v] = svd(w3+b2(numEigs*2+1:numEigs*3,:));    p3 = u * eye(numEigs,length(icIdx)) * v'; 
                
                b1 = b1 + A*w;
                b2 = b2 + (w - [p1;p2;p3]);
                ro1 = min(ro1 * 2,20/parcelNum * max(k/2,1));
                ro2 = min(ro2 * 4,40/(numEigs *3) * max(k/2,1));
                
                c1 = w1'*w1 - eye(length(icIdx)); 
                c2 = w2'*w2 - eye(length(icIdx)); 
                c3 = w3'*w3 - eye(length(icIdx));
                c4 = sum(sum(abs(A*w)))/parcelNum;
                fprintf('%d:%f: c1:%f c2:%f c3:%f c4:%f\n',iter,k, sum(sum(abs(c1))),sum(sum(abs(c2))),sum(sum(abs(c3))),c4);

                iter = iter + 1;
        end
    
        for i = 1:length(icIdx)
            filename = sprintf('./results_multi/IC_%d_1i_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X1'*w1(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename);
            filename = sprintf('./results_multi/IC_%d_2i_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X2'*w2(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename);
            filename = sprintf('./results_multi/IC_%d_3i_%f_%s.nii.gz',icIdx(i),k,subjects_012{subjectIdx});
            s = X3'*w3(:,i);
            new_ic = mask;
            new_ic.img(mask.img>0) = s;
            save_nii(new_ic,filename); 
        end
    end
    end    
    fprintf('\n %d finishing \n',subjectIdx);
end

