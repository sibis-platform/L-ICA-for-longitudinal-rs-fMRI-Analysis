function generateSlopes(icIdx,outputDir)

addpath('../Seed2Vox/nifti/');

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


mask = load_nii('/fs/neurosci01/djk/ncanda/atlas/sri24/restingstate/baseline/melodic/results_6mm_d25_dx/mask.nii.gz');

%% fslmerge
cmd =  sprintf('fslmerge -t ./slopes/slope_merged/slope_%d_0.05.nii.gz',icIdx);
cmd1 = sprintf('fslmerge -t ctrl_%d_0.05.nii.gz',icIdx);
cmd2 = sprintf('fslmerge -t etoh_%d_0.05.nii.gz',icIdx);
for subjectIdx = 1:length(subjects_012)
    group = strcmp(data_baseline_selected(strcmp(data_baseline_selected.subject,subjects_012{subjectIdx}),:).dx,'ctrl');
    if (group == 1)
        cmd1 = sprintf('%s ./%s/slope_%d_0.050000_%s', cmd1, outputDir, icIdx, subjects_012{subjectIdx});
    else
        cmd2 = sprintf('%s ./%s/slope_%d_0.050000_%s', cmd2, outputDir, icIdx, subjects_012{subjectIdx});
    end
    
    cmd = sprintf('%s ./%s/slope_%d_0.050000_%s', cmd, outputDir, icIdx, subjects_012{subjectIdx});

end

system(cmd);
% system(cmd1);
% system(cmd2);

end
