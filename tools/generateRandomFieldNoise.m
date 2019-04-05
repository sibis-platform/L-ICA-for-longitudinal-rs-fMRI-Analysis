function noiseData = generateRandomFieldNoise(t,mask,maskVoxNum)
%     system('rm noise*.nii');
%     system('rm all.nii');
%     system('rm 3dFWHMx*');
%     
%     cmd = sprintf('3dClustSim -nxyz %d %d %d -dxyz 2 2 2 -fwhm 17.28 -ssave:blurred noise.nii -iter %d -mask /fs/neurosci01/djk/ncanda/atlas/sri24/restingstate/baseline/melodic/results_6mm_d25_dx/mask.nii.gz',x,y,z,t);
%     system(cmd);
%     
%     system('fslmerge -t all.nii.gz noise*.nii');

    noiseData = zeros(t,maskVoxNum);    
    for i = 1:t
        idx = randi(1000) - 1;
        filename = sprintf('noise/noise_%6.6d.nii',idx);
        noise = load_nii(filename);
        noiseData(i,:) = noise.img(mask.img > 0);
    end
end