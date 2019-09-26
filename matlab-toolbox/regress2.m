function [B,BINT]=regress2(Y,X,alpha)
%        [B,BINT]=regress2(Y,X,[alpha])
% Linear regression of Y on X in which the confidence interval accounts for
% autocorrelation in record Y by using an AR(1) process as the null
% hypothesis rather than white noise as in the standard regress() function.
% Uses common method for estimating effective degrees of freedom in
% autocorrelated data (Bartlett, J. Roy. Stat. Soc. 98, 536-543, 1935;
% Mitchell et al., Climatic Change, WMO Technical Note 79, 1966; Bretherton
% et al., J. Climate 12, 1990-2009, 1999). Requires corr() and tinv()
% functions from Matlab Statistics Toolbox. Does not reproduce all of the
% functionality of regress(). See "help regress".
% Example [y is AR(1) process]:
% r=0.8; x=(1:100)'; y=zeros(100,1); u=randn(99,1); for n=2:100, y(n)=r*y(n-1)+u(n-1); end
% [b bint]   =  regress(y,[ones(size(x)) x]) % fit to y=b(1)+b(2)*x
% [b2 bint2] = regress2(y,[ones(size(x)) x])
% Ian Eisenman, 2010

if  nargin < 2
    error('regress2:input','regress2(Y,X,[alpha]): requires 2-3 input arguments.');
elseif nargin == 2
    alpha=0.05; % confidence level is 100*(1-alpha)
end
if ~isvector(Y) || size(Y,2)>1
    error('regress2:input','regress2(Y,X): Y must be a column vector.');
end
if size(Y,1)~=size(X,1)
    error('regress2:input','regress2(Y,X): X and Y must have same number of rows.')
end

% Code based on Matlab's regress().
[n,ncolX] = size(X);
% Use QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);
p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
B = zeros(ncolX,1);
B(perm) = R \ (Q'*Y);

if nargout >= 2
    % Find a confidence interval for each component of x
    r=corr(Y(1:end-1),Y(2:end));    % Estimate of AR(1) parameter using lag-1 autocorrelation coefficient
    neff=n*(1-r)/(1+r);             % Effective degrees of freedom (Mitchell et al 1966)
    nu = neff-p;                    % Residual degrees of freedom
    yhat = X*B;                     % Predicted responses at each data point.
    r = Y-yhat;                     % Residuals.
    rmse = norm(r)/sqrt(nu);        % Root mean square error.
    tval = tinv((1-alpha/2),nu);
    se = zeros(ncolX,1);
    se(perm,:) = rmse*sqrt(sum(abs(R\eye(p)).^2,2));
    BINT = [B-tval*se, B+tval*se];
end
