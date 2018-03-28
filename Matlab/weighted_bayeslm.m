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
%X2 = bsxfun(@minus,X,mean(X,1));
%if do_weights
 %   Y = Y-nanmean(Y,1);
  %  Y = Y./nanstd(Y,0,1);
%end
wY = Y;
[n,p] = size(X); % number of observations and parameters
%stdX = std(X,1) * sqrt(n);
%meanX = mean(X,1);
%In = ~~eye(n);
%X(~isfinite(X)) = 1;

%w = ones(n,1);
%lambda_i = gamrnd(1,1);
%ScoresL = Scores.*lambda_i;

%% Different feature types get different lambdas
%featureIDs(iCtrl) = {'Intercept'};
bind(1,:)  = strcmp(featureIDs,'Intercept');
bind(2,:)  = strcmp(featureIDs,'Group');
bind(3,:)  = strcmp(featureIDs,'Feature');
bind(4,:)  = strcmp(featureIDs,'Donor');
bind(5,:) = strcmp(featureIDs,'Feature:Group');
bind(6,:) = strcmp(featureIDs,'Feature:Donor');
bind(7,:)  = strcmp(featureIDs,'Subject');
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
%limit = lambda_i./(1/2+exp(nF/lambda_i^2)) + 0.5;
%X(:,~bind(1,:)) = bsxfun(@minus,X(:,~bind(1,:)),meanX(:,~bind(1,:)));
%X(:,~bind(1,:)) = bsxfun(@rdivide,X(:,~bind(1,:)),stdX(:,~bind(1,:)));
wX = X;

%% Initialize variables
rng('default');
tau_vector = rand(1,p);
D_tau_squared = diag(tau_vector);
XtX = wX'*wX;
%norm_betaG = zeros(2,p);
%pXtX = 1./XtX;
%Cov = diag(pXtX)';
%[norm_betaG,~,ic] = unique(round(Cov,2),'stable');
%nK = numel(norm_betaG);
%pXtX = diag(Cov);
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
XXtYMR = X(Ymissing(~MNR),:)*X(Ymissing(~MNR),:)';

%% Lambdas
lambda_ridge = rand(1,p);
%nu = rand(1,p);
lambda_lasso = rand(1,p);
%xi = rand(1,p);
%warning('off');
%% Gibbs sampler
for i = 1:iNumIter
    %index = 1:iNumIter;
    %sigma = sqrt(sigma2);
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
    beta_estimate = mvnrnd(A\wX'*wY,C)+1e-20;
    %beta_estimate(bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:)) = mvnrnd(A(bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:),bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:))\wX(:,bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:))'*wY,C(bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:),bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:)));%sampler(wX./sigma, A./sigma2, wY./sigma, C)';
    %if nFG > 0 || nFD > 0
     %   beta_estimate(bind(5,:)|bind(6,:)) = mvnrnd(A(bind(5,:)|bind(6,:),bind(5,:)|bind(6,:))\wX(:,bind(5,:)|bind(6,:))'*wY,C(bind(5,:)|bind(6,:),bind(5,:)|bind(6,:)));
    %end
    %% Sample sigma_squared from inverse gamma distribution
	sigma2_shape = (n-1+p)/2;%(n-1+p)/2 + 250/(1+nF)^2; %Play with this prior
    residuals = Y - (beta_estimate*X')';
	sigma2_scale = (residuals'*residuals./2 + (lambda_lasso(:).*beta_estimate(:))'/D_tau_squared*beta_estimate'./2 + lambda_ridge(:)'.*beta_estimate(:)'*beta_estimate'./2);
	sigma2 = 1./gamrnd(sigma2_shape,1./sigma2_scale+0.01);
    
    %% Sample 1/tau^2 from GIG distribution
    %for iii = 1:7
     %   norm_betaG(1,bind(iii,:)) = norm(beta_estimate(1,bind(iii,:)));
      %  tt = lambda_lasso(1,bind(iii,:));
       % if ~isempty(tt); norm_betaG(2,bind(iii,:)) = tt; end
    %end
    tau_shape = sqrt((lambda_lasso.^2.*sigma2)./beta_estimate.^2);%sqrt((norm_betaG(2,norm_betaG(1,:)~=0).^2.*sigma2)./norm_betaG(1,norm_betaG(1,:)~=0).^2);
    %tau_shape(isinf(tau_shape)) = 1e10;
    tau_scale = lambda_lasso.^2;%norm_betaG(2,norm_betaG(1,:)~=0).^2;
    tau_vector(1,:) =  gigrnd(-0.5, tau_shape, tau_scale);
    %t = zeros(1,p);
    %t(1,norm_betaG(1,:)~=0) = gigrnd(-0.5, tau_shape, tau_scale);      
    %for iii = 1:7
     %   tau_vector(1,bind(iii,:)) = t(1,bind(iii,:));
        
        %% Sample grouped lambda_ridge from gamma distributions
        %lambda_ridge(1,bind(iii,:)) = gamrnd(1, nALL(iii)/2/sigma2.*beta_estimate(1,bind(iii,:)).^2);
    %end
    
    %% Sample lambdas from gamma distributions
    %lambda_lasso(1,:) = sqrt(gamrnd(1+nG+nF+nD+nS, sum(1./tau_vector(1,bind(1,:)|bind(2,:)|bind(3,:)|bind(4,:)|bind(7,:)),2)+0.5));
    %lambda_lasso(1,:) = sqrt(gamrnd(p, sum(1./tau_vector,2)/2+0.5));
    %lambda_lasso(1,:) = sqrt(1./exprnd(1+tau_vector));
    lambda_lasso(1,:) = sqrt(gamrnd(p,sum(tau_vector,2).*2+1));
    %lambda_lasso(1,bind(1,:)) = 1e-3;
    %lambda_ridge(1,:) = gamrnd(p/2, 1/2/sigma2.*sum(beta_estimate.^2,2)+0.5);
    %lambda_ridge(1,:) = gamrnd((1+nF+nG+nD+nS)/2, 1/2/sigma2.*sum(beta_estimate.^2,2)+2);
    lambda_ridge(1,:) = gamrnd(1,1./(beta_estimate.^2./2./sigma2+3));
    %nu = 1 ./ exprnd(1+1./lambda_ridge);
    %lambda_lasso(1,bind(1,:)) = 1e-3;
    %xi = 1./ exprnd(1+1./lambda_lasso);
    %if nFG > 0 || nFD > 0
        %lambda_lasso(1,bind(5,:)|bind(6,:)) = sqrt(gamrnd(p, sum(tau_vector,2)/2+1));
     %   lambda_ridge(1,bind(5,:)|bind(6,:)) = gamrnd(p/2, 1/sum(beta_estimate.^2,2)/2/sigma2+1);
      %  lambda_ridge(1,bind(6,:)) = gamrnd(p, 1/sigma2.*sum(beta_estimate.^2,2)+1);
    %end
    %lambda_ridge(bind(1,:)) = gamrnd(1+nD*100, 1/2/sigma2.*sum(beta_estimate.^2,2));
    %lambda_ridge(bind(2,:)|bind(3,:)|bind(4,:)) = gamrnd((nD+1)/2, 1/2/sigma2.*sum(beta_estimate.^2,2));
    %lambda_ridge(bind(6,:)) = gamrnd(p/2, 1/2/sigma2.*sum(beta_estimate.^2,2) + 1/2/sigma2.*sum(beta_estimate(bind(1,:)|bind(6,:)).^2,2) + 1/2/sigma2.*sum(sum(beta_estimate(bind(1,:)|bind(3,:))'.^2*beta_estimate(bind(1,:)|bind(4,:)).^2,2),1));
    %lambda_ridge(bind(5,:)) = gamrnd(p/2, 1/2/sigma2.*sum(beta_estimate.^2,2));
%{
    %% Weak heridity implemented via tau and lambda_ridge for interaction coefficients
    if nFG > 0
        %tau_vector(1,bind(5,:)) = t(5);
        tt = 1./repmat(tau_vector(bind(2,:))',1,nF+1)./2 + 1./repmat(tau_vector(bind(1,:)|bind(3,:)),nG,1)./2;
        tt(iCtrl,:) = [];
        tau_vector(bind(5,:)) = tau_vector(bind(5,:)) .* 1./tt(:)';
        
        tt = repmat(lambda_ridge(bind(2,:))',1,nF+1)./2 + repmat(lambda_ridge(bind(1,:)|bind(3,:)),nG,1)./2;
        tt(iCtrl,:) = [];
        lambda_ridge(bind(5,:)) = lambda_ridge(bind(5,:)) .* tt(:)';
    end
    if nFD > 0
        %tau_vector(1,bind(6,:)) = t(6);
        tt = 1./repmat(tau_vector(bind(4,:))',1,nF+1)./2 + 1./repmat(tau_vector(bind(1,:)|bind(3,:)),nD,1)./2;
        tau_vector(bind(6,:)) = tau_vector(bind(6,:)) .* 1./tt(:)';
        
        tt = repmat(lambda_ridge(bind(4,:))',1,nF+1)./2 + repmat(lambda_ridge(bind(1,:)|bind(3,:)),nD,1)./2;
        lambda_ridge(bind(6,:)) = lambda_ridge(bind(6,:)) .* tt(:)';
    end
%}
    D_tau_squared = diag(tau_vector);
        
    if i > iBurn
        beta_posterior(ii,:) = beta_estimate;
        sigma_squared(ii) = sigma2;
        impY(:,ii) = Y(Ymissing);
        ii = ii + 1;
    end
    if do_weights
        %% Add weights to outlier observations based on Mascot Scores
        %lambda_i = 1.5+gamrnd(p,(2/sigma2)^2+0.001); %1.5+gamrnd(0.5+nF/20,2000/nF/sigma2+0.1);
        %ScoresL = Scores.*lambda_i^2;
        %limit = 1.345;%1/(lambda_i./(1/2+exp(lambda_i^2/nF)) + 0.5);
        %rTest = abs(residuals) > limit;%Decrease this??
        %w = ones(n,1);
        %w(rTest) = ScoresL(rTest);
        %S = wX*C*wX'*oo/n*p/3;
        %for iii = 1:n; r(iii) = 1/(1+residuals(iii)^2/2/sigma2 + S(iii,iii)/2/sigma2); end
        r = 1./(0.5+residuals.^2./2./sigma2);%+S./2./sigma2);% + sum(S,2)./2./sigma2);
        s = binornd(1,Scores);
        %r = Scores./(residuals.^2+S(In))./2./sigma2;
        %w = 1+gamrnd(2.*Scores, 1+r);
        w = 1+s+gamrnd(s+0.5,r);
        %w = 3.*Scores+3.*Scores./ max(1, abs(r));%(1+gamrnd(Scores+Scores.*abs(residuals./nanstd(residuals,0,1)),1/2/sigma2+0.01)).^2;
        wY = w.*Y;
        wX(:,~bind(1,:)) = bsxfun(@times,X(:,~bind(1,:)),w);
        XtX = wX'*wX;
        %Cov = 1./diag(XtX)';
    end
end

%% Imputed values from imputed distributions
impY = nanmean(impY,2);
Y(Ymissing) = impY;

%% Stats
%beta_posterior(:,~bind(1,:)) = bsxfun(@times,beta_posterior(:,~bind(1,:)),stdX(1,~bind(1,:)));
SEMs = nanstd(beta_posterior,0,1);
beta_estimate = nanmean(beta_posterior,1);
intercept = beta_estimate(bind(1,:));
yfit = (beta_estimate*X')';
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
%warning('on');
end