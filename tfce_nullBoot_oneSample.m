function [pValue,tfceStat] = tfce_nullBoot_oneSample(Y,h0,UseParFor)
% [pValue,tfceStat] = tfce_nullBoot_oneSample(Y,[h0],[UseParFor])
% Returns TFCE statistics and associated pValues for a one-sample...
% ... design testing the null hypothesis that observations are sampled...
% ... from a distribution with a mean of h0 (set to zero by default). ...
% ... The null distribution is computed by a bootstrapping procedure...
% ... that calculates the TFCE statistic after independently flipping...
% ... the sign of each observation at random (p=0.5).
% By Sam Berens (s.berens@sussex.ac.uk)
%
% INPUTS:
%    Y: Observed data, a cell array of NIfTI file names or a 4D numeric...
%       ... NOTE: if a numeric is provided, the last (4th) dimension...
%       ... must delineate independent observations.
%    h0: [Optional] Mean of the sampling distribution under the null.
%    UseParFor: [Optional] Parallelises bootstrapping using a parfor loop.
% OUTPUTS:
%    pValue: A 3D numeric of pValues.
%    tfceStat: A 3D numeric of TFCE statistics.

%% Check inputs
if nargin < 1
    error('tfce_nullBoot_oneSample calls require at least 1 input.');
end
if (nargin < 2) || isempty(h0)
    h0 = 0;
end
if nargin < 3
    UseParFor = true;
end

%% Hyperparamters
nBoot = 1e4;

%% Determine the number of observations (n)...
% ... and convert Y from a cell to a mat if needed.
if iscell(Y)
    if size(Y,1)==1
        Y = Y';
    end
    n = numel(Y);
    V = spm_vol(Y);
    Y = cellfun(@(vv)spm_read_vols(vv),V,'UniformOutput',false);
    Y = cell2mat(permute(Y,[2,3,4,1]));
else
    n = size(Y,4);
end

%% Subtract h0 from Y
Y = Y - h0;

%% Set X and H
X = ones(n,1);
H = 1;

%% Compute the tfce statistic
tfceStat = tfce_Xcon(Y,X,H);

%% Vectorise the tfce statistic
mask = tfceStat > 0;
tfceStat_vec = tfceStat(mask);

%% Run bootstrapping
nullStat = nan([numel(tfceStat_vec),nBoot]);
fprintf('Bootstrapping...%c',10)
if UseParFor
    nullStat = ParBoot(n,X,Y,H,mask,nullStat);
else
    nullStat = SerBoot(n,X,Y,H,mask,nullStat);
end

%% Compute the pValues
pValue_vec = mean(repmat(tfceStat_vec,1,nBoot)>nullStat,2);
pValue = nan(size(tfceStat));
pValue(mask) = pValue_vec;

return

function [nullStat] = ParBoot(n,X,Y,H,mask,nullStat)
nBoot = size(nullStat,2);
parfor iBoot = 1:nBoot
    cY = reshape(randsample([-1,1],n,true),[1,1,1,n]) .* Y;
    nullStat_vol = tfce_Xcon(cY,X,H);
    nullStat(:,iBoot) = nullStat_vol(mask);
end
return

function [nullStat] = SerBoot(n,X,Y,H,mask,nullStat)
nBoot = size(nullStat,2);
for iBoot = 1:nBoot
    cY = reshape(randsample([-1,1],n,true),[1,1,1,n]) .* Y;
    nullStat_vol = tfce_Xcon(cY,X,H);
    nullStat(:,iBoot) = nullStat_vol(mask);
    if mod(iBoot,79)==0
        fprintf('... %06.2f%% complete%c',(iBoot/nBoot)*100,10)
    end
end
return