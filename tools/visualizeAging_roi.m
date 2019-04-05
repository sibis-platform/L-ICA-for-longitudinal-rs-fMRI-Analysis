clear;

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

mask  = load_nii('./masks/sri24_functional_parcellation_100.nii.gz');
roiIdx = 27;
icIdx = 2;
regNum = 0.05;
M = max(mask.img(:));

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

S = zeros(M,sum(mask.img(:) > 0));
p = mask.img(mask.img > 0);
xc = XC(mask.img > 0);
yc = YC(mask.img > 0);
zc = ZC(mask.img > 0);

for i = 1:max(mask.img(:))
	idx = (p == i);
	S(i,idx) = 1;
	mxc = mean(xc(idx));
	myc = mean(yc(idx));
	mzc = mean(zc(idx));
    
	dis = (xc - mxc).^2 + (yc - myc).^2 + (zc - mzc).^2;
    
	S(i,:) = exp(-dis / 50);
	S(i,:) = S(i,:) / norm(S(i,:));
end   

x = [];
y = [];

beta_ctrl = [];
beta_etoh = [];
for i = 1:size(subjects_012,1)
    subjectName = subjects_012{i};
	filename1 = sprintf('./results_multi/IC_%d_1_%f_%s.nii.gz',icIdx,regNum,subjectName);
    filename2 = sprintf('./results_multi/IC_%d_2_%f_%s.nii.gz',icIdx,regNum,subjectName);
    filename3 = sprintf('./results_multi/IC_%d_3_%f_%s.nii.gz',icIdx,regNum,subjectName);
    
    if (exist(filename1) == 0) || (exist(filename2) == 0) || (exist(filename3) == 0)
        continue;
    end
    group = strcmp(data_baseline_selected(strcmp(data_baseline_selected.subject,subjects_012{i}),:).dx,'ctrl');
    
    
    img1 = load_nii(filename1);
	img2 = load_nii(filename2);
	img3 = load_nii(filename3);    
    
    age1 = table2array(data_baseline(strcmp(data_baseline.subject,subjects_012{i}),6));
	age2 = table2array(data_f1y(strcmp(data_f1y.subject,subjects_012{i}),6));
	age3 = table2array(data_f2y(strcmp(data_f2y.subject,subjects_012{i}),6));
    
%     measure1 = mean(img1.img(mask.img == roiIdx));
%     measure2 = mean(img2.img(mask.img == roiIdx));
% 	  measure3 = mean(img3.img(mask.img == roiIdx));
    measure1 = S(roiIdx,:) * img1.img(mask.img > 0 ) ;
    measure2 = S(roiIdx,:) * img2.img(mask.img > 0 ) ;
    measure3 = S(roiIdx,:) * img3.img(mask.img > 0 ) ;   
    
    x = [x;age1;age2;age3];
    y = [y;measure1;measure2;measure3];

    if (group > 0.5)
        plot([age1,age2,age3],[measure1,measure2,measure3],'lineWidth',1,'color',[1,0.7,0.7]);
        beta_ctrl = [beta_ctrl, [1,age1;1,age2;1,age3]\[measure1;measure2;measure3]];
    else
        plot([age1,age2,age3],[measure1,measure2,measure3],'lineWidth',2,'color',[0.7,0.7,1]);
        beta_etoh = [beta_etoh, [1,age1;1,age2;1,age3]\[measure1;measure2;measure3]];
    end
    
    hold on;
        
    fprintf('Finish %d\n',i);
end

beta_ctrl = mean(beta_ctrl,2);
beta_etoh = mean(beta_etoh,2);

plot([12;24],[1,12;1,24]*beta_ctrl,'r','lineWidth',3);
plot([12;24],[1,12;1,24]*beta_etoh,'b','lineWidth',3);