clear;

data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline_ex0.csv');
data_baseline_selected = data_baseline((strcmp(data_baseline.dx,'ctrl') | strcmp(data_baseline.dx,'heavy') | ...
                                        strcmp(data_baseline.dx,'moderate'))&...
                                       data_baseline.b_restingstate == 1, :);
                                   
data_f1y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_1y_ex0.csv');
data_f1y_selected = data_f1y((strcmp(data_f1y.dx,'ctrl') | strcmp(data_f1y.dx,'heavy') |...
                                          strcmp(data_f1y.dx,'moderate'))&...
                                          data_f1y.b_restingstate == 1, :);
                                      
                                      
data_f2y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_2y_ex0.csv');
data_f2y_selected = data_f2y((strcmp(data_f2y.dx,'ctrl') | strcmp(data_f2y.dx,'heavy') |...
                                          strcmp(data_f2y.dx,'moderate'))&...
                                          data_f2y.b_restingstate == 1, :);

subjects_02 = intersect(data_baseline_selected.subject,data_f2y_selected.subject);
subjects_01 = intersect(data_baseline_selected.subject,data_f1y_selected.subject);
subjects_012 = intersect(data_f2y_selected.subject,subjects_01);

data = data_baseline_selected(ismember(data_baseline_selected.subject,subjects_012),:);

mask = load_nii('/fs/neurosci01/djk/ncanda/atlas/sri24/restingstate/baseline/melodic/results_6mm_d25_dx/mask.nii.gz');
%parcellation = load_nii('masks/sri24_functional_parcellation_100.nii.gz');
parcellation = load_nii('../melodic_group_ICA/25IC/regions_simplified/clusters_index18_simplified.nii.gz');
M = max(parcellation.img(:));   

ctrl_num = 0;
etoh_num = 0;
icIdx = 18;
regNum = 0.01;

filename = sprintf('/fs/neurosci01/qingyuz/rsfmri/melodic_group_ICA/25IC/stats/thresh_zstat%d.nii.gz',icIdx);
foreground_mask = load_nii(filename);
foreground_mask = foreground_mask.img(mask.img > 0);
foreground_mask = (foreground_mask > 0.01);
XC = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
YC = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
ZC = zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));

for i = 1:size(mask.img,1)
        for j = 1:size(mask.img,2)
            for k = 1:size(mask.img,3)
                XC(i,j,k) = i;
                YC(i,j,k) = j;
                ZC(i,j,k) = k;
            end
        end
end

%S = zeros(M,sum(mask.img > 0));
p = parcellation.img(mask.img > 0);
xc = XC(mask.img > 0);
yc = YC(mask.img > 0);
zc = ZC(mask.img > 0);

sidx = 0;
for i = 1:max(p)
    idx = (p == i);

    dice = sum(idx & foreground_mask) / sum(idx);
    if (dice < 0.5)
        continue;
    end
    
    sidx = sidx + 1;
    S(sidx,:) = idx;
	S(sidx,:) = S(sidx,:) / sum(S(sidx,:));
    
end   
M = sidx;

response = zeros(length(subjects_012),M);
for subjectIdx = 1:length(subjects_012)
    filename1 = sprintf('./results_multiv3/IC_%d_1_%f_%s.nii.gz',icIdx,regNum,subjects_012{subjectIdx});
	filename2 = sprintf('./results_multiv3/IC_%d_2_%f_%s.nii.gz',icIdx,regNum,subjects_012{subjectIdx});
	filename3 = sprintf('./results_multiv3/IC_%d_3_%f_%s.nii.gz',icIdx,regNum,subjects_012{subjectIdx});
                    
	if (exist(filename1) == 0) || (exist(filename2) == 0) || (exist(filename3) == 0)
        continue;
    end  
        
    raw1 = load_nii(filename1);
    raw2 = load_nii(filename2);    
    raw3 = load_nii(filename3);
    
	age1 = table2array(data_baseline(strcmp(data_baseline.subject,subjects_012{subjectIdx}),6));
	age2 = table2array(data_f1y(strcmp(data_f1y.subject,subjects_012{subjectIdx}),6));
	age3 = table2array(data_f2y(strcmp(data_f2y.subject,subjects_012{subjectIdx}),6));
    
    img1 = S*raw1.img(mask.img > 0 );
    img2 = S*raw2.img(mask.img > 0 );
    img3 = S*raw3.img(mask.img > 0 );
    
    X = [eye(M),diag(ones(1,M)*age1);...
         eye(M),diag(ones(1,M)*age2);...
         eye(M),diag(ones(1,M)*age3)];
     
    beta = X\[img1;img2;img3];

    response(subjectIdx,:) = beta(M+1:2*M);
    subjectIdx
end

filename1 = sprintf('response%d.mat',icIdx);
save(filename1, 'response');

%%
filename1 = sprintf('response%d.mat',icIdx);
load(filename1);

score = readtable('score.csv');
design = [strcmp(data.sex,'M') - strcmp(data.sex,'F')];
design = [design,strcmp(data.scanner,'ge') - strcmp(data.scanner,'siemens')];
design = [design,data.visit_age];
design = [design,data.visit_age .* (strcmp(data.sex,'M') - strcmp(data.sex,'F'))];
design = [design,score.Var2];

tbl = design;

for j = 1:M
    if size(response,1) == 1
        [b,dev,stats{j}] = glmfit(tbl,response{j});
    else
        [b,dev,stats{j}] = glmfit(tbl,response(:,j));
    end
end


filename2 = sprintf('tbl%d.mat',icIdx);
filename3 = sprintf('stats%d.mat',icIdx);

save(filename2, 'tbl');
save(filename3, 'stats');

%% glm fdr

p = zeros(1,M);
for i = 1:M
    p(i) = stats{i}.p(end);
end

[p,idx] = sort(p);
fdr_alpha = 0.1 / M;

for i = 1:M
    if p(i) <= fdr_alpha * i
        idx(i)
    end
end

min(p)
0.1/M


