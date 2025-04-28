function [pValue,tfceStat] = tfce_nullBoot_corr(Y,x,UseParFor)
% [pValue,tfceStat] = tfce_nullBoot_corr(Y,x,[UseParFor])
% Returns TFCE statistics and associated pValues for a correlational...
% ... design testing the null hypothesis that the covariance between...
% ... the input vector x and the observed data Y is zero. ...
% ... The null distribution is computed by a bootstrapping procedure...
% ... that calculates the TFCE statistic after independently flipping...
% ... each observation (Y) around the mean at random (p=0.5).
% By Sam Berens (s.berens@sussex.ac.uk)
%
% INPUTS:
%    Y: Observed data, a cell array of NIfTI file names or a 4D numeric...
%       ... NOTE: if a numeric is provided, the last (4th) dimension...
%       ... must delineate independent observations.
%    x: An [n,1] predictor variable, where n is the number of observations.
%    UseParFor: [Optional] Parallelises bootstrapping using a parfor loop.
% OUTPUTS:
%    pValue: A 3D numeric of pValues.
%    tfceStat: A 3D numeric of TFCE statistics.

%% Check inputs
if nargin < 2
    error('tfce_nullBoot_corr calls require at least 2 inputs.');
end
if iscell(Y)
    if size(Y,1)==1
        Y = Y';
    end
    if numel(x) ~= numel(Y)
        error(['Input x should be an [n,1] vector, ',...
            'where n is the number of observations in Y.']);
    end
else
    if numel(x) ~= size(Y,4)
        error(['Input x should be an [n,1] vector, ',...
            'where n is the number of observations in Y.']);
    end
end
if nargin < 3
    UseParFor = true;
end

%% Hyperparamters
nBoot = 1e4;

%% Determine the number of observations (n)...
% ... and convert Y from a cell to a mat if needed.
if iscell(Y)
    n = numel(Y);
    V = spm_vol(Y);
    Y = cellfun(@(vv)spm_read_vols(vv),V,'UniformOutput',false);
    Y = cell2mat(permute(Y,[2,3,4,1]));
else
    n = size(Y,4);
end

%% Set X and H ...
% ... Ensure X is mean centred and has a unit length
if size(x,1) ~= n
    x = x';
end
x = x - mean(x);
x = x ./ norm(x);
X = x;
H = 1;

%% Ensure Y is mean centred and has a unit length
Y = Y - mean(Y,4);
Y = num2cell(Y,4);
Y = cellfun(@(v) v./norm(squeeze(v)),Y,'UniformOutput',false);
Y = cell2mat(Y);

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