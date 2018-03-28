%% Imputation function. 
%Assumes that missing values are missing due to low
%abundance. Imputes missing values (NaNs) from a normal distribution with a
%mean and std derived from the minimum feature values with a centre shifted
%by -1.6 * stdev.
function [Y] = Impute(X, N)
Y = X;

switch N
    case 'none'   %Only takes features that don't have missing values
        for i = 1:size(X,1)
            if any(isnan(X(i,:)))
                Y(i,:) = NaN;
            end
        end
        
    case 'MI'
        vals = nanmean(X,2);
        sigmas = nanstd(X,0,2);
        mus = vals - (1.6 .* sigmas);
        for i = 1:size(X,1)
            ImputedVals = (0.3 * sigmas(i)) .*...
                randn(size(X(i,isnan(X(i,:))),2),1) + mus(i);
            Y(i,isnan(Y(i,:))) = ImputedVals;
        end
        
    otherwise
        figure('Name', 'Total and Imputed Distributions')
        for i = 1:size(X,2)
            vals = nanmin(X,[],2);
            sigma = nanstd(vals,0,1);
            mu = nanmean(vals,1) - (1.6 * sigma);
            ImputedVals = (0.3 * sigma) .*...
                randn(size(X(isnan(X(:,i)),i),1),1) + mu;
            Y(isnan(Y(:,i)),i) = ImputedVals;
    
            %Set subplot numbers before running program
            subplot(floor(sqrt(size(X,2)))+1,floor(sqrt(size(X,2)))+1,i)
            histogram(Y(:,i),30)
            hold on;
            histogram(ImputedVals,20)
            hold off;
            axis([-10 40 0 20000]);
        end
end

end