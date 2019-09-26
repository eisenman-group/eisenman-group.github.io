function [RHO,PVAL]=corr2(X,Y)
%        [RHO,PVAL]=corr2(X,Y)
% Linear correlation (RHO) between X and Y in which the p-value (PVAL)
% accounts for autocorrelation in both records by using AR(1) processes as
% the null hypothesis rather than white noise as in the standard corr()
% function. Uses common method for estimating effective degrees of freedom
% in autocorrelated data (Bartlett, J. Roy. Stat. Soc. 98, 536-543, 1935;
% Mitchell et al., Climatic Change, WMO Technical Note 79, 1966; Bretherton
% et al., J. Climate 12, 1990-2009, 1999). Requires corr() and tcdf()
% functions from Matlab Statistics Toolbox. Does not reproduce all of the
% functionality of corr(). See "help corr".
% Example [x and y are AR(1) processes]:
% rx=0.7; x=zeros(100,1); u=randn(99,1); for n=2:100, x(n)=rx*x(n-1)+u(n-1); end
% ry=0.9; y=zeros(100,1); u=randn(99,1); for n=2:100, y(n)=ry*y(n-1)+u(n-1); end
% [r p]=corr(x,y)
% [r2 p2]=corr2(x,y)
% Ian Eisenman, 2010

if  nargin < 2
    error('corr2:input','corr2(X,Y): requires two input arguments');
end
if ~isvector(Y)
    error('corr2:input','corr2(X,Y): Y must be a vector.');
end
if ~isvector(X)
    error('corr2:input','corr2(X,Y): X must be a vector.');
end
if sum(size(Y)~=size(X))
    error('corr2:input','regress2(Y,X): X and Y must be same size.');
end

dx=X-mean(X);
dy=Y-mean(Y);

RHO = (dx./norm(dx))' * (dy./norm(dy)); % correlation coefficient

if nargout==2
    rx=corr(X(1:end-1),X(2:end));          % estimate of AR(1) parameter using lag-1 autocorrelation coefficient
    ry=corr(Y(1:end-1),Y(2:end));          % estimate of AR(1) parameter using lag-1 autocorrelation coefficient
    neff=length(Y)*(1-rx*ry)/(1+rx*ry);    % effective degrees of freedom (Mitchell et al 1966)
    t=abs(RHO)*sqrt( (neff-2)/(1-RHO^2) ); % T-value
    PVAL=2*tcdf( -t,neff-2 );              % Student's T cumulative distribution function, two-tailed test
end
