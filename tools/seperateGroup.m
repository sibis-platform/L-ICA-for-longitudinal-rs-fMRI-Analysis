clear;
addpath('/fs/neurosci01/qingyuz/rsfmri/Seed2Vox/nifti/');

data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline_ex0.csv');
data_baseline_selected = data_baseline((strcmp(data_baseline.dx,'ctrl') | strcmp(data_baseline.dx,'ctrl') | ...
                                        strcmp(data_baseline.dx,'ctrl'))&...
                                       data_baseline.b_restingstate == 1, :);
                                   
data_f1y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_1y_ex0.csv');
data_f1y_selected = data_f1y((strcmp(data_f1y.dx,'ctrl') | strcmp(data_f1y.dx,'moderate') |...
                                          strcmp(data_f1y.dx,'heavy'))&...
                                          data_f1y.b_restingstate == 1, :);
                                      
                                      
data_f2y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_2y_ex0.csv');
data_f2y_selected = data_f2y((strcmp(data_f2y.dx,'ctrl') | strcmp(data_f2y.dx,'moderate') |...
                                          strcmp(data_f2y.dx,'heavy'))&...
                                          data_f2y.b_restingstate == 1, :);

subjects_02 = intersect(data_baseline_selected.subject,data_f2y_selected.subject);
subjects_01 = intersect(data_baseline_selected.subject,data_f1y_selected.subject);
subjects_012 = intersect(data_f2y_selected.subject,subjects_01);

%%
data = data_baseline_selected(ismember(data_baseline_selected.subject,subjects_012),:);
N = size(data,1);

flag = 1;
p_age_max = 0;
p_sex_max = 0;
p_scanner_max = 0;

while flag < 1000
    idx = randperm(N);
    
    data1 = data(idx(1:N/2),:);
    data2 = data(idx(1+N/2:end),:);
    
    %%
    ob = [sum(strcmp(data1.sex,'M')), sum(strcmp(data1.sex,'F'));...
          sum(strcmp(data2.sex,'M')), sum(strcmp(data2.sex,'F'))];

    bn1 = sum(strcmp(data.sex,'M')) / N; 
    bn0 = sum(strcmp(data.sex,'F') )/ N;

    ex = [size(data1,1) * bn1, size(data1,1) * bn0;...
          size(data2,1) * bn1, size(data2,1) * bn0];
  
    chis = (ob - ex).^2 ./ ex;
    p_sex = 1 - chi2cdf(sum(chis(:)),1);
    
    %%
    [h_age,p_age] = ttest2(data1.visit_age,data2.visit_age);
    
    %%
    
    ob = [sum(strcmp(data1.scanner,'ge')), sum(strcmp(data1.scanner,'siemens'));...
          sum(strcmp(data2.scanner,'ge')), sum(strcmp(data2.scanner,'siemens'))];

    bn1 = sum(strcmp(data.scanner,'ge')) / N; 
    bn0 = sum(strcmp(data.scanner,'siemens') )/ N;

    ex = [size(data1,1) * bn1, size(data1,1) * bn0;...
          size(data2,1) * bn1, size(data2,1) * bn0];
  
    chis = (ob - ex).^2 ./ ex;
    p_scanner = 1 - chi2cdf(sum(chis(:)),1);
    
    if (p_scanner > 0.2) && (p_age > 0.2) && (p_sex > 0.2) && ...
       (p_scanner > p_scanner_max && p_age > p_age_max && p_sex > p_sex_max)    
        p_scanner_max = max(p_scanner_max,p_scanner);
        p_age_max = max(p_age_max,p_age);
        p_sex_max = max(p_sex_max,p_sex);
        best_idx = idx;
    end
    %% 
    flag = flag + 1;
end
