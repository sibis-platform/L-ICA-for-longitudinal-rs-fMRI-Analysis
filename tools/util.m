clear;
addpath('../Seed2Vox/nifti');

path = '../vox_statistical_test/palm_results_ctrl_0.200000_100.000000_onesample';
icIdx = [1:25];

cmd = sprintf('rm %s/finding*',path);
system(cmd);

threshold = 0.0025;
clusterSize = 0;

for i = icIdx
    filename = sprintf('%s/%d_tfce_tstat_fwep_c1.nii',path,i);
    if exist(filename)
        p = load_nii(filename);
    
        p = p.img(p.img > 0);
        %min(p)
        if sum(p <= threshold)>clusterSize
            fprintf('bingo: %d_c1 %d %f !\n',i,sum(p <= threshold),min(p));
            p = load_nii(filename);
            location = (p.img <= threshold & p.img > 0);
            p.img(~location) = 0;
            p.img(location) = 1;
            p.img = flip(p.img,1);
            filename = sprintf('%s/finding_ctrl_%d_tfce_c1.nii',path,i);
            save_nii(p,filename);
        end
    end
    
    filename = sprintf('%s/%d_tfce_tstat_fwep_c2.nii',path,i);
    if exist(filename)
        p = load_nii(filename);
    
        p = p.img(p.img > 0);
        %min(p)
        if sum(p <= threshold)>clusterSize
            fprintf('bingo: %d_c2 %d %f !\n',i,sum(p <= threshold),min(p));
            p = load_nii(filename);
            location = (p.img <= threshold & p.img > 0);
            p.img(~location) = 0;
            p.img(location) = 1;
            p.img = flip(p.img,1);
            filename = sprintf('%s/finding_ctrl_%d_tfce_c2.nii',path,i);
            save_nii(p,filename);
        end
    end
    
%     filename = sprintf('%s/%d_tfce_tstat_fwep_c3.nii',path,i);
%     if exist(filename)
%         p = load_nii(filename);
%     
%         p = p.img(p.img > 0);
%         %min(p)
%         if sum(p <= threshold)>clusterSize
%             fprintf('bingo: %d_c3 %d !\n',i,sum(p <= threshold));
%             p = load_nii(filename);
%             location = (p.img <= threshold & p.img > 0);
%             p.img(~location) = 0;
%             p.img(location) = 1;
%             p.img = flip(p.img,1);
%             filename = sprintf('%s/finding_ctrl_%d_tfce_c3.nii',path,i);
%             save_nii(p,filename);
%         end
%     end
%     
%         filename = sprintf('%s/%d_tfce_tstat_fwep_c4.nii',path,i);
%     if exist(filename)
%         p = load_nii(filename);
%     
%         p = p.img(p.img > 0);
%         %min(p)
%         if sum(p <= threshold)>clusterSize
%             fprintf('bingo: %d_c4 %d !\n',i,sum(p <= threshold));
%             p = load_nii(filename);
%             location = (p.img <= threshold & p.img > 0);
%             p.img(~location) = 0;
%             p.img(location) = 1;
%             p.img = flip(p.img,1);
%             filename = sprintf('%s/finding_ctrl_%d_tfce_c4.nii',path,i);
%             save_nii(p,filename);
%         end
%     end
end
