%% Implementation of the Devroye (2014) algorithm for sampling from 
%% the generalized inverse Gaussian (GIG) distribution
%
%% Modified by Venkatesh Mallikarjun to allow P vector input on 19/07/2017
%
% function X = gigrnd(p, a, b)
%
% The generalized inverse Gaussian (GIG) distribution is a continuous
% probability distribution with probability density function:
%
% p(x | p,a,b) = (a/b)^(p/2)/2/besselk(p,sqrt(a*b))*x^(p-1)*exp(-(a*x + b/x)/2)
%
% Parameters:
%   p \in Real, a > 0, b > 0
%
% See Wikipedia page for properties:
% https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution 
%
% This is an implementation of the Devroye (2014) algorithm for GIG sampling.
% 
% Returns:
%   X      = random variates [size(p,2)] from the GIG(p, a, b)
%  
% References:
% L. Devroye
% Random variate generation for the generalized inverse Gaussian distribution 
% Statistics and Computing, Vol. 24, pp. 239-246, 2014.
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2015
function X = gigrnd(P, a, b)

if isscalar(P)
    P = repmat(P,size(a));
end
%% Setup -- we sample from the two parameter version of the GIG(alpha,omega)
lambda = P;
omega = sqrt(a.*b);
alpha = sqrt(omega.^2 + lambda.^2) - lambda;

%% Find t
x = -psi(1, alpha, lambda);
t = NaN(size(lambda));
t((x >= 0.5) & (x <= 2)) = 1;
t(x > 2) = sqrt(2 ./ (alpha(x > 2) + lambda(x > 2)));
t(x < 0.5) = log(4./(alpha(x < 0.5) + 2.*lambda(x < 0.5)));
%{
if((x >= 0.5) && (x <= 2))
    t = 1;
elseif(x > 2)
    t = sqrt(2 / (alpha + lambda));
elseif(x < 0.5)
    t = log(4/(alpha + 2.*lambda));     
end
%}
%% Find s
x = -psi(-1, alpha, lambda);
s = NaN(size(lambda,2));
s((x >= 0.5) & (x <= 2)) = 1;
s(x > 2) = sqrt(4./(alpha(x > 2).*cosh(1) + lambda(x > 2)));
s(x < 0.5) = min(1./lambda(x < 0.5), log(1 + 1./alpha(x < 0.5) + sqrt(1./alpha(x < 0.5).^2+2./alpha(x < 0.5))));
s = s(:);
%{
if((x >= 0.5) && (x <= 2))
    s = 1;
elseif(x > 2)
    s = sqrt(4./(alpha.*cosh(1) + lambda));
elseif(x < 0.5)
    s = min(1./lambda, log(1 + 1./alpha + sqrt(1./alpha.^2+2./alpha)));
end
%}
%% Generation
eta = -psi(t, alpha, lambda);
zeta = -dpsi(t, alpha, lambda);
theta = -psi(-s, alpha, lambda);
xi = dpsi(-s, alpha, lambda);
p = 1./xi;
p(~isfinite(p)) = 1;
r = 1./zeta;
r(~isfinite(r)) = 1;
td = t - r.*eta;
sd = s - p.*theta;
q = td + sd;

sampleSize = size(p,2);
X = zeros(1,sampleSize);
for sample = 1:sampleSize
    done = false;
    iter = 1;
    while(~done)
        U = rand(1); 
        V = rand(1); 
        W = rand(1);
        if(U < (q(sample) / (p(sample) + q(sample) + r(sample))))
            X(sample) = -sd(sample) + q(sample).*V;
        elseif(U < ((q(sample) + r(sample)) / (p(sample) + q(sample) + r(sample))))
            X(sample) = td(sample) - r(sample).*log(V);
        else
            X(sample) = -sd(sample) + p(sample).*log(V);
        end

        %% Are we done?
        f1 = exp(-eta - zeta.*(X(sample)-t(sample)));
        f2 = exp(-theta + xi.*(X(sample)+s(sample)));
        if((W*g(X(sample), sd(sample), td(sample), f1(sample), f2(sample))) <= exp(psi(X(sample), alpha(sample), lambda(sample))))
            done = true;
        elseif iter > 10000
                error('gigrnd could not converge');
        else
            iter = iter + 1;
        end
    end
end

%% Transform X back to the three parameter GIG(p,a,b)
X = exp(X) .* (lambda ./ omega + sqrt(1 + (lambda./omega).^2));
X = X ./ sqrt(a./b);

end

function f = psi(x, alpha, lambda)
    f = -alpha.*(cosh(x(:)) - 1) - lambda.*(exp(x(:)) - x(:) - 1);
end

function f = dpsi(x, alpha, lambda)
    f = -alpha.*sinh(x(:)) - lambda.*(exp(x(:)) - 1);
end

function f = g(x, sd, td, f1, f2)

a = 0;
b = 0;
c = 0;
if((x >= -sd) && (x <= td))
    a = 1;
elseif(x > td)
    b = f1;
elseif(x < -sd)
    c = f2;   
end;

f = a + b + c;

end