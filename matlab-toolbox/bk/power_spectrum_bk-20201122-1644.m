function [P,s,ci]=power_spectrum(y,smooth)

%    [P,s,ci]=power_spectrum(y,[smooth])
%
% Usage: Creates power spectral estimate P(s), where s is freqency in 
%   units of (y data spacing)^-1.
% Outputs: P = power; s = frequency [cycles/dt]; ci = 0.95 confidence interval.
%   The confidence interval is based on white noise. ci given at all s if
%   the multitaper method is used (a 2xN array); otherwise it is a 2x1 array.
% Inputs:
%  If smooth=0, periodogram (absolute value of Fourier transform) is computed.
%  If smooth>O, estimated power spectral density is computed by averaging over 
%   1+2*smooth frequency bands (default: smooth=2). Note that spectral density times
%   N=length(y) is equivalent to a smoothing of the periodogram, so smooth=0 and 
%   smooth>0 results differ by a factor of length(y). The number of degrees of 
%   freedom is 2+4*smooth.
%  If smooth<0, a multi-taper spectrum is computed with pmtm using nw=abs(smooth)+1 
%   for the number of windows ("time bandwidth product") (typically smooth=-2); 
%   this gives approximately 2+4*smooth degrees of freedom.
%  Note: Spectrum has max period of record length and min (Nyquist) period of twice 
%   the sampling interval.
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
% Ian Eisenman, 2006


% Longer example 1: compare periodogram, frequency band averaged power
%   spectral density estimate, and multitaper method
if 0 % not to be evaluated in power_spectrum: for pasting into command window
    % create time series
    T=100; dt=0.01; t=dt:dt:T; T1=3; T2=7; N=length(t);
    y=cos(t*(2*pi/T1))+cos(t*(2*pi/T2))+2*randn(1,N);
    % create spectra by 3 methods
    [P1,s01,ci1]=power_spectrum(y,0); s1=s01/diff(t(1:2)); % s01=cycles/dt; s1=cycles/t
    [P2,s02,ci2]=power_spectrum(y,2); s2=s02/diff(t(1:2));
    [P3,s03,ci3]=power_spectrum(y,-2); s3=s03/diff(t(1:2));
    % create plots
    subplot(2,1,1) % plot time series
    plot(t,y,t,cos(t*(2*pi/T1)),t,cos(t*(2*pi/T2))), axis tight
    subplot(2,3,4) % plot 3 spectral density estimates
    loglog(s1,P1/N,s2,P2,'g',s3,P3,'r')
    ylabel('power density'), xlabel('frequency (cycles/t)')
    y0=10^2.5; x1=10^-1.8; x2=x1*2; x3=x1*4; % plot fixed confidence intervals
    hold on, plot([x1 x1],y0*ci1,'b',x1,y0,'bo')
    plot([x2 x2],y0*ci2,'g',x2,y0,'go')
    plot([x3 x3],y0*mean(ci3),'r',x3,y0,'ro'), hold off
    set(gca,'ylim',10.^[-2 3.5],'xlim',10.^[-2 1.5])
    subplot(2,3,5) % plot spectra vs period instead of frequency
    plot(1./s1,P1/N,'b.-',1./s2,P2,'g.-',1./s3(2:end),P3(2:end),'r.-')
    axis tight, set(gca,'yscale','log','xlim',[0 15],'ylim',10.^[-2 3.5]); grid on
    xlabel('period (t)'), legend('periodogram','spectral density','multi-taper',4)
    subplot(2,3,6) % plot multitaper spectrum with shaded confidence intervals
    fill([s3(2:end); flipud(s3(2:end))],[ci3(2:end,1).*P3(2:end);flipud(ci3(2:end,2).*P3(2:end))],...
        0.8*[1 1 1],'edgecolor','none'); hold on, plot(s3,P3,'r'), hold off
    set(gca,'xscale','log','yscale','log','ylim',10.^[-2 3.5],'xlim',10.^[-2 1.5])
    title('multitaper spectrum')
    disp(['multitaper confidence interval limits vary with s by '...
            num2str(std(ci3)./mean(ci3)*100,'%0.3g ') '%'])
end
% end of longer example

% Longer example 2: compute percent variance explained
if 0
    % create time series
    t=0:0.01:200;
    y1=8.3*cos(t*2*pi/0.5);
    y2=3.1*cos(t*2*pi/20+1.3);
    % create and plot power spectrum (should have peak at s=1/period)
    [P,s,ci]=power_spectrum(y1+y2+randn(size(y1)),0);
    % find s grid interfaces; note that s is evenly spaced
    si=[s-diff(s(1:2))/2; s(end)+diff(s(1:2))/2];
    subplot(1,2,1)
    loglog(s/diff(t(1:2)),P), axis tight
    ylabel('power density'), xlabel('frequency (cycles/t)')
    subplot(1,2,2)
    plot(s/diff(t(1:2)),P), axis tight, set(gca,'xscale','log')
    % percent of variance in y1
    P2=P.*diff(si);
    disp([var(y2)/(var(y1+y2)) 3.1^2/(8.3^2+3.1^2) ...
        sum(P2(s/diff(t(1:2))<1))/sum(P2)])
end
% 


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
    P=abs(Y(2:N/2+1))'.^2; % Y(N/2+1) is value at Nyqist frequency
    s=(1:N/2)'/N;
    
elseif smooth>0 % spectral density estimate
    if smooth~=floor(smooth), disp('Error: smooth must be integer if smooth>0'), return; end

    % from Wunsch 12.864 course at MIT, Spring 2006
    % power density spectral estimate
    A=fft(y)/length(y);
    N=2*floor(length(y)/2);
    dt=1;
    T=N*dt;
    s=(1:N/2)'/T;

%    B=floor(nu/4);
    B=smooth; % B=(nu-2)/4, where nu is number of degrees of freedom
    PD=zeros(1,N/2+1)*nan;
    for n=(1+B):(N/2+1-B)
        PD(n) = N/(2*B+1) * sum(abs(A(n-B:n+B)).^2); % (Wunsch notes 11.3)
    end
    P=PD(2:end)';
    
else % smooth<0: multi-taper spectral density estimate
%    [P,s,ci]=pmtmPH(y,1,abs(smooth)+1,0,length(y));
    [P0,ci0,s]=pmtm(y,abs(smooth)+1,length(y),1);
    P=P0/2; ci=ci0./repmat(P0,1,2); % normalization and meaning of ci conventions
    P=P(2:end); ci=ci(2:end,:); s=s(2:end); % discard frequency=s=0 value
end

% remove NAN values at edges
g=find(isnan(P)); P(g)=[]; s(g)=[]; if smooth<0, ci0(g,:)=[]; end
