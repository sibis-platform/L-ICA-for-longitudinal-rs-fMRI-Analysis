function [whitenedData,w_init,ref,numSamples] = whiteningData(mask,raw,IC,icIdx,numEigs)

    %% re-organize data and normalization
    maskVoxNum = sum(mask.img(:) > 0);
    rawData = zeros(size(raw.img,4),maskVoxNum);
    regression_x = zeros(25,maskVoxNum);

    for t = 1:size(raw.img,4)
        raw_img = squeeze(raw.img(:,:,:,t));
        rawData(t,:) = raw_img(mask.img > 0)';
    end
    for ic = 1:25
        ic_img = squeeze(IC.img(:,:,:,ic));
        regression_x(ic,:) = ic_img(mask.img > 0)';
    end

    normalizedData = rawData;

    %% pca and whitening
    numSamples = maskVoxNum;
    pcaData = normalizedData - repmat(mean(normalizedData,2),1,size(normalizedData,2));
    [whitenedData] = fastica(pcaData,'only','white', 'lastEig', numEigs, 'verbose','off');

    %% dual regression initialization
    %w_init = glmfit(whitenedData',regression_x(icIdx,:),'normal','constant','off');
    w_init = whitenedData' \ regression_x(icIdx,:)';
    scale = norm(w_init);
    w_init = w_init/scale;
    ref = regression_x(icIdx,:)/scale;
end