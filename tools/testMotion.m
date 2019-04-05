clear;
clf;

data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline_ex0.csv');
data_baseline_selected = data_baseline((strcmp(data_baseline.dx,'ctrl') | strcmp(data_baseline.dx,'ctrl') | ...
                                        strcmp(data_baseline.dx,'ctrl'))&...
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

%%
score = readtable('../palm/score.csv');
score = score(score.Var2 < 80,:);
subjects_012 = intersect(subjects_012, score.subjects_012);

%%
motion = readtable('motion.csv');

A = [1,0;1,1;1,2];

all = [];
for i = 1:length(subjects_012)
    visit = motion.Var3(strcmp(subjects_012{i},motion.Var1));
    mp = motion.Var4(strcmp(subjects_012{i},motion.Var1));
    
    all = [all;mp(1:3)];
    beta = A \ mp(1:3);
    ms(i) = beta(2);
end

save ms.mat ms
data = data_baseline_selected(ismember(data_baseline_selected.subject,subjects_012),:);

[h,p_sex] = ttest2(ms(strcmp(data.sex,'M')) , ms(strcmp(data.sex,'F'))) ;

[b,dev,stat] = glmfit(data.visit_age,ms);
[h,p_onesample] = ttest(ms);