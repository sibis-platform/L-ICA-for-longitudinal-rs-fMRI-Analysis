function wrapperLICA(subjectIdx) 
    if ~isdeployed
        addpath('/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi_25/FastICA_25/');
        addpath('/fs/neurosci01/qingyuz/rsfmri/Seed2Vox/nifti/');
        addpath('/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi_25/tools/');
    end
    if isdeployed
        [compThreads, count] = sscanf(getenv('NSLOTS'), '%d');
        if count == 1
            fprintf('NSLOTS=%d\n', compThreads);
            warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
            autoCompThreads = maxNumCompThreads(compThreads);
        end
    end

    mask = load_nii('/fs/neurosci01/qingyuz/rsfmri/dual_regression_3y_ctrl/results_8mm_d25_dx/mask.nii.gz');
    IC = load_nii('/fs/neurosci01/qingyuz/rsfmri/dual_regression_3y_ctrl/results_8mm_d25_dx/melodic_IC.nii.gz');
    

    
    filename = sprintf('./filenames/%s.mat',subjectIdx);
    load(filename);
    filename = sprintf('./age/age_%s.mat',subjectIdx);
    load(filename);

    if (length(filenames) ~= length(age))
        fprintf('%s data not match\n',subjectIdx);
        return;
    end
    
    outputPath = './test';
    Options.icIdx = [1:11,13:22,24,25];
    Options.ICAOption = 1;
    Options.k1 = 0.02;
    Options.k2 = 0.2;
    Options.parcellation = '/fs/neurosci01/qingyuz/rsfmri/Longitudinal_ICA_multi_25/masks/sri24_functional_parcellation_100.nii.gz';
    
    singleSubjectMultiICA_general(subjectIdx, filenames, mask, IC, age, outputPath, Options);
end
