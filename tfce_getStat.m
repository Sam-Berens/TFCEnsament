function [tfceStat] = tfce_getStat(M)
% [tfceStat] = tfce_getStat(M)
% Returns TFCE statistics for an image of test statistics.
% By Sam Berens (s.berens@sussex.ac.uk)
% 
% M: A 3D numeric of test statistics (e.g., t- or F-values).
% tfceStat: A 3D numeric of TFCE statistics.

%% Hyperparamters
dh = 0.1;
hPower = 2;
ePower = 0.5;

%% Determin the range of hights that we must loop over
maxh = max(M,[],'all');
h = (0:dh:maxh)';
nSteps = numel(h);

%% Loop through each hight to calculate the extent
E = zeros([size(M),nSteps]);
sizeE = size(E);
for iStep = 1:nSteps
    B = M > h(iStep);
    CC = bwconncomp(B,26);
    k = cellfun(@numel,CC.PixelIdxList)';
    nClusters = numel(k);
    for iC = 1:nClusters
        [i1,i2,i3] = ind2sub(CC.ImageSize,CC.PixelIdxList{iC});
        IL = sub2ind(sizeE,i1,i2,i3,ones(size(i1)).*iStep);
        E(IL) = k(iC);
    end
end

%% Loop through each hight to calculate the tfce stat
tfceStat = zeros(size(M));
for iStep = 1:nSteps
    tfceStat = tfceStat + (E(:,:,:,iStep).^ePower).*(h(iStep)^hPower).*dh;
end

return