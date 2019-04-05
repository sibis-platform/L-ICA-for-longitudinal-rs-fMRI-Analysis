function [whitenedData,w_init,ref,numSamples] = whiteningDataMulti(mask,rawData,regression_x,icIdx,numEigs)

    %% re-organize data and normalization
    maskVoxNum = sum(mask.img(:) > 0);
    normalizedData = rawData;

    %% pca and whitening
    numSamples = maskVoxNum;
    pcaData = normalizedData - repmat(mean(normalizedData,2),1,size(normalizedData,2));
    [whitenedData] = fastica(pcaData,'only','white', 'lastEig', numEigs, 'verbose','off');

    %% dual regression initialization
    w_init = zeros(numEigs,length(icIdx));
    ref = zeros(length(icIdx),maskVoxNum);
    for i = 1:length(icIdx)
        %w_init(:,i) = glmfit(whitenedData',regression_x(icIdx(i),:),'normal','constant','off');
        w_init(:,i) = whitenedData' \ regression_x(icIdx(i),:)';
        scale = norm(w_init(:,i));
        w_init(:,i) = w_init(:,i)/scale;
        ref(i,:) = regression_x(icIdx(i),:)/scale;
    end
end