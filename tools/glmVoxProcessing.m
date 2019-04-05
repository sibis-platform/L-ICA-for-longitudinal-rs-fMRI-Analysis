function glmVoxProcessing(icIdx,fileDir)

filename = sprintf('./%s/slope_%d_0.05.nii.gz',fileDir,icIdx);
if (exist(filename))
    cmd = sprintf('fslchfiletype NIFTI %s',filename);
    system(cmd);
end
return

% addpath('../Seed2Vox/nifti/');
% addpath('../Seed2Vox/palm-alpha109/');
% 
% data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline_ex0.csv');
% data_baseline_selected = data_baseline((strcmp(data_baseline.dx,'ctrl') | ...
%                                         strcmp(data_baseline.dx,'heavy'))&...
%                                        data_baseline.b_restingstate == 1, :);
%                                    
% data_f1y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_1y_ex0.csv');
% data_f1y_selected = data_f1y((strcmp(data_f1y.dx,'ctrl') | ...
%                                           strcmp(data_f1y.dx,'heavy'))&...
%                                           data_f1y.b_restingstate == 1, :);
%                                       
%                                       
% data_f2y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_2y_ex0.csv');
% data_f2y_selected = data_f2y((strcmp(data_f2y.dx,'ctrl') | ...
%                                           strcmp(data_f2y.dx,'heavy'))&...
%                                           data_f2y.b_restingstate == 1, :);
% 
% subjects_02 = intersect(data_baseline_selected.subject,data_f2y_selected.subject);
% subjects_01 = intersect(data_baseline_selected.subject,data_f1y_selected.subject);
% subjects_012 = intersect(data_f2y_selected.subject,subjects_01);
% 
% data = data_baseline_selected(ismember(data_baseline_selected.subject,subjects_012),:);
% 
% N = size(data,1);
% regNum = 0.05;
% 
% %%
% design = [strcmp(data.sex,'M') - strcmp(data.sex,'F')];
% design = [design,strcmp(data.scanner,'ge') - strcmp(data.scanner,'siemens')];
% design = [design,data.visit_age];
% 
% design = [design,(data.race == 1) - (data.race == 6)];
% design = [design,(data.race == 2) - (data.race == 6)];
% design = [design,(data.race == 3) - (data.race == 6)];
% design = [design,(data.race == 4) - (data.race == 6)];
% design = [design,(data.race == 5) - (data.race == 6)];
% design = [design,strcmp(data.dx,'ctrl') - strcmp(data.dx,'heavy')];
% 
% %%
% 
% filename = sprintf('./%s/slope_%d_0.05.nii',fileDir,icIdx);
% img = load_nii(filename);
% mask = load_nii('./masks/mask.nii');
% 
% result_t = mask;
% result_p = mask;
% 
% for i = 1:size(img.img,1)
%     for j = 1:size(img.img,2)
%         for k = 1:size(img.img,3)
%             if (mask.img(i,j,k) > 0)
%                 response = squeeze(img.img(i,j,k,:));
%                 [b,dev,stats] = glmfit(design,response);
%                 result_t.img(i,j,k) = stats.t(10);
%                 result_p.img(i,j,k) = stats.p(10);
%             end
%         end
%         i
%         j
%     end
% end
% 
% filenamet = sprintf('%d_glm_t.nii.gz',icIdx);
% filenamep = sprintf('%d_glm_p.nii.gz',icIdx);
% save_nii(result_t,filenamet);
% save_nii(result_p,filenamep);

end

