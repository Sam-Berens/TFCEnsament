function [pValue,tfceStat] = tfce_nullBoot_permuteRows(Y,X,H)
% [pValue,tfceStat] = tfce_nullBoot_permuteRows(Y,X,H)
% Returns TFCE statistics and associated pValues for a GLM contrast (H),...
% ... by permuting rows of the design matrix (X) to bootstrap the null...
% ... distribution.
% By Sam Berens (s.berens@sussex.ac.uk)
% 
% Y: Observed data, as a cell array of NIfTI file names or a 4D numeric.
% X: GLM design matrix.
% H: GLM contrast matrix.
% pValue: A 3D numeric of pValues.
% tfceStat: A 3D numeric of TFCE statistics.

%% Hyperparamters
nBoot = 1e4;

%% Determine the number of observations (n)
n = size(X,1);

%% Compute the tfce statistic
tfceStat = tfce_Xcon(Y,X,H);

%% Vectorise the tfce statistic
mask = tfceStat > 0;
tfceStat_vec = tfceStat(mask);

%% Run bootstrapping
nullStat = nan([numel(tfceStat_vec),nBoot]);
fh = waitbar(0,'Bootstrapping...');
for iBoot = 1:nBoot
    p = randperm(n)';
    cX = X(p,:);
    nullStat_vol = tfce_Xcon(Y,cX,H);
    nullStat(:,iBoot) = nullStat_vol(mask);
    if mod(iBoot,17)==0
        waitbar(iBoot/nBoot,fh);
    end
end
close(fh);

%% Compute the pValues
pValue_vec = mean(repmat(tfceStat_vec,1,nBoot)>nullStat,2);
pValue = nan(size(tfceStat));
pValue(mask) = pValue_vec;

return