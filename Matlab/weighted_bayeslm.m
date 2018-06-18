function [mdl] = weighted_bayeslm(X,Y,featureIDs,do_weights,Scores,MNR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - weighted_bayeslm
% Author: Venkatesh Mallikarjun
%
% Description: 
%   Function to perform Bayesian Elastic Net with parameter-specific 
%   weights calculated using standard deviations, and weighted 
%   residuals if required.
%
% Input:  
%       X - matrix of predictors
%       Y - vector of responses
%       featureIDs - cell row vector of strings denoting parameter types.
%       do_weights - Determines whether to do residual weighting.
%       Scores - Vector of Mascot scores or other indicator of peptide ID
%       confidence.
%       MNR - n x 1 binary vector denoting which observations in Y are
%       missing non-randomly
%
% Output: 
%       1) Matrix of inferred regressors (mdl.B)
%       2) average estimate of regressors (mdl.beta_estimate)
%       3) model intercept (mdl.intercept)
%       4) Mean Squared Error of model fit (mdl.MSE)
%       5) Standard Error of regressor estimates (mdl.SEPred)
%       6) t-scores for each beta estimate (mdl.tscores)
%       7) p-values for beta_estimate (mdl.Pvalues)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wY = Y;
[n,p] = size(X); % number of observations and parameters

bind(1,:)  = strcmp(featureIDs,'Intercept');
bind(2,:)  = strcmp(featureIDs,'Group');
bind(3,:)  = strcmp(featureIDs,'Feature');
bind(4,:)  = strcmp(featureIDs,'Donor');
bind(5,:) = strcmp(featureIDs,'Feature:Group');
bind(6,:) = strcmp(featureIDs,'Feature:Donor');
bind(7,:)  = strcmp(featureIDs,'Subject');

% Scale features (0<= X <= 1) for when ContGroup is used
minX = min(X(:,~bind(1,:)),[],1);
X(:,~bind(1,:)) = bsxfun(@minus,X(:,~bind(1,:)),minX);
maxX = max(X(:,~bind(1,:)),[],1);
X(:,~bind(1,:)) = bsxfun(@rdivide,X(:,~bind(1,:)),maxX);
%nU = sum(any(bind,2),1);
%nK2 = sum(nK,1);
nG = sum(bind(2,:),2);
nF = sum(bind(3,:),2);
nD = sum(bind(4,:),2);
%nFG = sum(bind(5,:),2);
%nFD = sum(bind(6,:),2);
nS = sum(bind(7,:),2);
%nALL = [1,nG,nF,nD,nFG,nFD,nS];
iNumIter = max(1000,25*max(nS,nG*(nD+1))+nF);
iBurn = ceil(iNumIter/2);
wX = X;

%% Initialize variables
rng('default');
tau_vector = rand(1,p);
D_tau_squared = diag(tau_vector);
XtX = wX'*wX;
w = ones(n,1);
beta_posterior = NaN(iNumIter-iBurn,p);
beta_estimate = randn(1,p);%(pXtX*(wX(~isnan(wY),:)'*wY(~isnan(wY))))'; %OLS fit
sigma_squared = NaN(iNumIter-iBurn, 1);
sigma2 = 1./gamrnd((n-1+p)/2,0.01);%nanmean((wY - wX*beta_estimate').^2,1)+1e-6; %OLS sigma_squared

%Imputation variables
impmin = nanmin(Y)-2;
Ymissing = find(isnan(wY));
nMissing = numel(Ymissing);
nMNR = sum(MNR,1);
nMR = sum(~MNR,1);
prop_MNR = nMNR/n;
if nMissing; alpha = prctile(Y,prop_MNR*100); end
impY = NaN(size(Ymissing,1),iNumIter-iBurn);
ii = 1;
XXtYMR = [ones(nMR,1),X(Ymissing(~MNR),:)]*[ones(nMR,1),X(Ymissing(~MNR),:)]';

%% Lambdas
lambda_ridge = rand(1,p);
lambda_lasso = rand(1,p);

%% Gibbs sampler
for i = 1:iNumIter
    %index = 1:iNumIter;
    %sigma = sqrt(sigma2);
    %b0 = nanmean(Y) - mean(beta_estimate*X');
    if nMissing
        if nMNR
            %Impute MNR missing values from truncated gaussian
            plo=normcdf(impmin/sigma2)+1e-100;
            phi=normcdf(alpha/sigma2)+1e-50;
            z=norminv(plo+(phi-plo).*rand(nMNR,1));
            wY(Ymissing(MNR)) = z*sigma2.*w(Ymissing(MNR));
        end
        if nMR
            %XXtYMR = wX(Ymissing(~MNR),:)*wX(Ymissing(~MNR),:)';
            B = sigma2./XXtYMR.*eye(nMR);
            wY(Ymissing(~MNR)) = mvnrnd(X(Ymissing(~MNR),:)*beta_estimate',B)'.*w(Ymissing(~MNR));%sampler(wX(Ymissing(~MNR),:)'./sigma,XXtYMR./sigma2,beta_estimate'./sigma,B);
        end
        Y(Ymissing) = wY(Ymissing)./w(Ymissing);
    end
    
    %% beta posterior (conditional)
    L = diag(lambda_ridge) + D_tau_squared;%sqrt(lambda_ridge'*lambda_ridge).*eye(p)+D_tau_squared;
    A = XtX + L;
    C = sigma2.*inv(A);
    beta_estimate = mvnrnd(A\wX'*wY,C);
    
    %% Sample sigma_squared from inverse gamma distribution
    sigma2_shape = (n-1+p)/2;%(n-1+p)/2 + 250/(1+nF)^2; %Play with this prior
    residuals = Y - (beta_estimate*X')';
    sigma2_scale = (residuals'*residuals./2 + (lambda_lasso(:).*beta_estimate(:))'/D_tau_squared*beta_estimate'./2 + lambda_ridge(:)'.*beta_estimate(:)'*beta_estimate'./2);
    sigma2 = 1./gamrnd(sigma2_shape,1./sigma2_scale+0.01);
    
    %% Sample 1/tau^2 from GIG distribution
    tau_vector(1,:) = 1./(1+gigrnd(0.5, lambda_lasso(:)./(4.*lambda_ridge(:).*sigma2), (lambda_ridge(:).*beta_estimate(:).^2)./sigma2)).^2;
    D_tau_squared = diag(tau_vector);
    
    %% Sample lambdas from gamma distributions
    lambda_lasso(1,:) = sqrt(gamrnd(p,sum(1./tau_vector,2)/2+1));
    lambda_ridge(1,:) = gamrnd(1,1./(beta_estimate.^2./2./sigma2+3));
        
    if i > iBurn
        beta_posterior(ii,:) = beta_estimate;
        sigma_squared(ii) = sigma2;
        impY(:,ii) = Y(Ymissing);
        ii = ii + 1;
    end
    if do_weights
        %% Add weights to outlier observations based on Mascot Scores
        r = 1./(0.5+residuals.^2./2./sigma2);%+S./2./sigma2);% + sum(S,2)./2./sigma2);
        s = binornd(1,Scores);
        w = 1+s+gamrnd(s+0.5,r);
        wY = w.*Y;
        wX(:,~bind(1,:)) = bsxfun(@times,X(:,~bind(1,:)),w);
        XtX = wX'*wX;
    end
end

%% Imputed values from imputed distributions
impY = nanmean(impY,2);
Y(Ymissing) = impY;

%% Stats
beta_posterior(:,~bind(1,:)) = bsxfun(@minus,beta_posterior(:,~bind(1,:)),minX);
beta_posterior(:,~bind(1,:)) = bsxfun(@rdivide,beta_posterior(:,~bind(1,:)),maxX);
SEMs = nanstd(beta_posterior,0,1);
beta_estimate = nanmean(beta_posterior,1);
yfit = (beta_estimate*X')';
intercept = mean(Y) - mean(yfit);
MSE = nanmean(sigma_squared);

%% return model
mdl.B = beta_posterior;
mdl.beta_estimate = beta_estimate;
mdl.intercept = intercept;
mdl.yfit = yfit;
mdl.residErr = MSE;
mdl.SEPred = SEMs;
mdl.Yimputed = Y;
mdl.tscores = abs(beta_estimate)./SEMs;
mdl.DoF = numel(beta_estimate(beta_estimate ~= 0));
mdl.Pvalues = (1 - tcdf(mdl.tscores', mdl.DoF - 1))' .* 2; %2-sided t-test
mdl.FeatureType = featureIDs;
end
