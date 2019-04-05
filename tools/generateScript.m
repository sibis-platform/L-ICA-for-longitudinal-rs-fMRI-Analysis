data_baseline = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_baseline.csv');
data_baseline_selected = data_baseline(strcmp(data_baseline.exceeds_bl_drinking, 'Y') &...
                                       data_baseline.R11 == 1 & data_baseline.b_structural == 1 & data_baseline.b_restingstate == 1, :);
                                   
data_f1y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_1y.csv');
data_f1y_selected = data_f1y(strcmp(data_f1y.exceeds_bl_drinking, 'Y') &...
                             data_f1y.R11 == 1 & data_f1y.b_structural == 1 & data_f1y.b_restingstate == 1, :);
                                     
data_f2y = readtable('/fs/neurosci01/djk/ncanda/group/scripts/design_matrix/ncanda_followup_2y.csv');
data_f2y_selected = data_f2y(strcmp(data_f2y.exceeds_bl_drinking, 'Y') &...
                             data_f2y.R11 == 1 & data_f2y.b_structural == 1 & data_f2y.b_restingstate == 1, :);
 
subjects_02 = intersect(data_baseline_selected.subject,data_f2y_selected.subject);
subjects_01 = intersect(data_baseline_selected.subject,data_f1y_selected.subject);
subjects_012 = intersect(data_f2y_selected.subject,subjects_01);

fid = fopen('run_ica_on_cluster_ex.sh','w');
sumNum = 1;
for subjectIdx = 1:length(subjects_012)
% 	out_filename_1 = sprintf('/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi/results_multi/IC_2_1_0.050000_%s.nii.gz',subjects_012{subjectIdx});
%     out_filename_2 = sprintf('/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi/results_multi/IC_2_2_0.050000_%s.nii.gz',subjects_012{subjectIdx});
%     out_filename_3 = sprintf('/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi/results_multi/IC_2_3_0.050000_%s.nii.gz',subjects_012{subjectIdx});
% 
% 	if (exist(out_filename_1)) && (exist(out_filename_2)) && (exist(out_filename_3))
%         fprintf('\n %d finishing \n',subjectIdx);
%         continue;
%     end
    
    filename1 = sprintf('/fs/ncanda-share/cases/%s/standard/baseline/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
	filename2 = sprintf('/fs/ncanda-share/cases/%s/standard/followup_1y/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
	filename3 = sprintf('/fs/ncanda-share/cases/%s/standard/followup_2y/restingstate/reslice/sri24_2mm/bold_noIntenCorr_4d_filtermotart_cleaned_bp_smooth6mm_stripped.nii.gz',...
                        subjects_012{subjectIdx});
	if (exist(filename1) == 0) || (exist(filename2) == 0) || (exist(filename3) == 0)
        continue;
    end
    
    fprintf(fid,'title=ica_%s\n',subjects_012{subjectIdx});
    fprintf(fid,'log=${title}.log\n');
    fprintf(fid,'cmd="/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi/ica.sh %s"\n',subjects_012{subjectIdx});
    fprintf(fid,'echo "${cmd}" | qsub -S /bin/bash -cwd -o ${log} -e ${log} -j y -pe smp 1 -l h_vmem=16G -N ${title}\n');
    
    sumNum = sumNum + 1;
end

fclose(fid);