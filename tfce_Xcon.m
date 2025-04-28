function [tfceStat,testStat,df1,df2] = tfce_Xcon(Y,X,H)
% [tfceStat,testStat,df1,df2] = tfce_Xcon(Y,X,H)
% Returns TFCE statistics for a GLM contrast.
% By Sam Berens (s.berens@sussex.ac.uk)
% 
% Y: Observed data, as a cell array of NifTi file names or a 4D numeric.
% X: GLM design matrix.
% H: GLM contrast matrix.
% tfceStat: A 3D numeric of TFCE statistics.
% testStat: A 3D numeric of t- or F-values (depending on the rank of H).
% df1: Contrast degrees of freedom (df1 = rank(H)).
% df2: Residual degrees of freedom (df2 = trace(M)).

%% Check the input Y to see whether it is a cell array of image names
% ... If so, load the images as a matrix
if iscell(Y)
    if size(Y)==1
        Y = Y';
    end
    V = spm_vol(Y);
    Y = cellfun(@(vv)spm_read_vols(vv),V,'UniformOutput',false);
    Y = cell2mat(permute(Y,[2,3,4,1]));
end

%% Ensure that the first dimension of Y relates to different observations
% ... and then reshape Y to collapse across spatial dimensions
if size(X,1)~=size(Y,1)
    Y = permute(Y,[4,1,2,3]);
end
volSize = size(Y,[2,3,4]);
nVox = prod(volSize);
Y = reshape(Y,[size(Y,1),nVox]);

%% OLS to compute the testStat 
B = pinv(X) * Y;
Bcov = pinv(X'*X);
P = X * pinv(X);
M = eye(size(P,1)) - P;
R = M * Y;
df1 = rank(H);
df2 = trace(M);
testStat = nan(nVox,1);
for iVox = 1:nVox
    mse = R(:,iVox)'*R(:,iVox) ./df2;
    C = Bcov.*mse;
    if df1 == 1
        testStat(iVox) = (H*B(:,iVox)) / sqrt(H*C*H');
    else
        testStat(iVox) = (H*B(:,iVox))'*(inv(H*C*H'))*(H*B(:,iVox)) / df1;
    end
end

%% Reshape the testStat and compute the tfce stat
testStat = reshape(testStat,volSize);
tfceStat = tfce_getStat(testStat);

return