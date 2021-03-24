function [P,s,ci]=power_spectrum(y,smooth)

%    [P,s,ci]=power_spectrum(y,[smooth])
%
% Usage: Creates power spectral density estimate P(s) for y(t), where s is freqency in
%   units of (dt)^-1 (where dt is the sampling period of the time series y.)
% Outputs: P = power; s = frequency [cycles/dt]; ci = 95% confidence interval.
%   The confidence interval is based on white noise. ci is given at all s if
%   the multitaper method is used (a 2xN array); otherwise it is a 2x1 array.
% Inputs:
%  If smooth=0, periodogram (proportional to absolute value of Fourier transform) is computed.
%  If smooth>O, frequency band averaged power spectral density estimate, computed by averaging over 
%   1+2*smooth frequency bands (default: smooth=2). Note that the spectral density 
%   estimate is equivalent to a smoothing of the periodogram. The number of degrees of 
%   freedom is 2+4*smooth.
%  If smooth<0, Thompson multitaper method power spectral density estimate, computed with pmtm using nw=abs(smooth)+1 
%   for the number of windows ("time bandwidth product") (typically smooth=-2); 
%   this gives approximately 2+4*smooth degrees of freedom.
%  Note: Spectrum has max period of record length and min (Nyquist) period of twice 
%   the sampling interval.
%  Normalization: The integral under the spectral density curve, sum(P)*diff(s(1:2)), 
%   is (approximately) equal to the signal variance, var(y).
%
% Example: Plot power spectrum with 95% confidence interval
%   % create time series
%   smooth=2; t=0:0.01:200; period=5;
%   y=cos(t*2*pi/period)+randn(size(t));
%   % create and plot power spectrum (should have peak at s=1/period)
%   [P,s,ci]=power_spectrum(y,smooth);
%   ci=mean(ci,1); % in case smooth<0 is used (does nothing when smooth>0)
%   loglog(s/diff(t(1:2)),P), axis tight
%   ylabel('power density'), xlabel('frequency (cycles/t)')
%   % add line indicating 95% confidence interval
%   xl=get(gca,'xlim'); yl=get(gca,'ylim');
%   x0=exp(log(xl(1))+diff(log(xl))*0.1); % 10% from left edge
%   y0=exp(log(yl(1))+diff(log(yl))*0.7); % 70% from bottom
%   hold on, plot([x0 x0],y0*ci,'k',x0,y0,'ko'), hold off
%
% Ian Eisenman, 2006, updated 2020

if nargin<1, disp('[P,s,ci]=power_spectrum(y,[smooth]); loglog(s/diff(t(1:2)),P)'), return; end
if nargin<2, smooth=2; end

if size(y,2)==1, y=y'; end
y=y-mean(y); % remove mean

% confidence interval for periodogram or power spectral density estimate
if smooth>=0
    % 0.95 conifidence interval for a chi-squared distribution,
    % cf., chi2confPH and pmtmPH.m
    %ci confidence interval
    %v  degrees of freedom
    % Approximation from Chamber's et al 1983; see Percival and Walden 1993 p.256
    v=2+4*smooth;
    ci=1./[(1-2./(9*v)+1.96*sqrt(2./(9*v))).^3 (1-2./(9*v)-1.96*sqrt(2./(9*v))).^3];
    % While this approximation is terrific for v>=6 (minimum smoothing), it
    % gives 2x too high an upper limit when v=2 (periodogram).
    % Manually correct this:
    if v==2, ci(2)=39.5; end
    % %Can use stats toolbox chi2inv for a more accurate approx; it's very
    % %similar to C83/PW93 result above
    % if exist('chi2inv')>0,  % need the stats toolbox for chi2inv
    %     cl=0.95; %cl confidence level
    %     ci=[v/chi2inv(1-(1-cl)/2,v) v/chi2inv((1-cl)/2,v)];
end

if smooth==0 % periodogram
    Y=fft(y);
    N=2*floor(length(y)/2); % number of elements in Y, made even if odd
    P=abs(Y(2:N/2+1))'.^2*2/N; % Y(N/2+1) is value at Nyqist frequency
    s=(1:N/2)'/N;
    
elseif smooth>0 % spectral density estimate
    if smooth~=floor(smooth), disp('Error: smooth must be integer if smooth>0'), return; end

    % from Carl Wunsch 12.864 course at MIT, Spring 2006
    % power density spectral estimate
    A=fft(y)/length(y);
    N=2*floor(length(y)/2);
    dt=1;
    T=N*dt;
    s=(1:N/2)'/T;
    B=smooth; % B=(nu-2)/4, where nu is number of degrees of freedom
    PD=zeros(1,N/2+1)*nan;
    for n=(1+B):(N/2+1-B)
        PD(n) = N/(2*B+1) * sum(abs(A(n-B:n+B)).^2); % (Wunsch notes 11.3)
    end
    P=2*PD(2:end)';
    
else % smooth<0: multi-taper spectral density estimate
    [P0,ci0,s]=pmtm(y,abs(smooth)+1,length(y),1);
    P=P0; ci=ci0./repmat(P0,1,2); % convention for meaning of ci
    P=P(2:end); ci=ci(2:end,:); s=s(2:end); % discard frequency=s=0 value
end

% remove NAN values at edges
g=find(isnan(P)); P(g)=[]; s(g)=[]; if smooth<0, ci0(g,:)=[]; end
