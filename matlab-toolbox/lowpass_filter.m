function y_filtered = lowpass_filter(y,freq,N)
% y_filtered = lowpass_filter(y,freq,[N])
% Lowpass filter
% y: unfiltered time series
% freq: cutoff frequency (cycles / time spacing of y)
% N: type of filter (default is N=10)
%   If N=0, uses simple fft method (fft, zero out frequencies outside range, then ifft)
%   If N>0, uses N-point butterworth filter
%   If N<0, uses N-order Chebyshev type I filter
% Example:
%  % 200 years sampled 10x per year, 3-year lowpass filter
%  y=rand(2001,1); t=0:1/10:200;
%  filter_N=10; lowpass_period=3;
%  yf=lowpass_filter(y,diff(t(1:2))/lowpass_period,filter_N);
%  subplot(2,1,1), plot(t,y,t,yf), xlabel('t'), ylabel('x')
%  [P,s,ci]=power_spectrum(y,3); Pf=power_spectrum(yf,3);
%  subplot(2,1,2), loglog(s/diff(t(1:2)),P,s/diff(t(1:2)),Pf), axis tight
%  hold on, loglog(1/lowpass_period*[1 1],[1e-10 1e10],'k--'), hold off, xlabel('f'), ylabel('P')
% Ian Eisenman, 2010

if nargin==0
    disp('y_filtered = lowpass_filter(y,freq,[N])')
end
if nargin==2
    N=10;
end

if 0 % longer example
    clf
    y=rand(2001,1); t=0:1/10:200; lowpass_period=3;
    filter_N=2; yf0=lowpass_filter(y,diff(t(1:2))/lowpass_period,filter_N);
    filter_N=6; yf1=lowpass_filter(y,diff(t(1:2))/lowpass_period,filter_N);
    filter_N=12; yf2=lowpass_filter(y,diff(t(1:2))/lowpass_period,filter_N);
    [P,s,ci]=power_spectrum(y,3); Pf0=power_spectrum(yf0,3); Pf1=power_spectrum(yf1,3); Pf2=power_spectrum(yf2,3);
    subplot(2,1,1), plot(t,y,t,yf0,t,yf1,'--',t,yf2,':'), xlabel('t'), ylabel('x')
    subplot(2,1,2), loglog(s/diff(t(1:2)),P,s/diff(t(1:2)),Pf0,s/diff(t(1:2)),Pf1,'--',s/diff(t(1:2)),Pf2,':'), axis tight
    hold on, loglog(1/lowpass_period*[1 1],[1e-10 1e10],'k--'), hold off, xlabel('f'), ylabel('P')
end

% remove mean
ym=mean(y);
y=y-ym;

if N==0 % fft method
    
    Y=fft(y);
    n=length(Y)-1; %first element of Y is sum(temp1)
    f=(1:n/2)/(n/2)*1/2; %half as long as Y
    ncf=length(find(f<freq));
    Y(ncf+1:n-ncf+1)=0; % remove high frequency variation
    y_filtered=real(ifft(Y)); %filtered data
    
elseif N>0 % butterworth method
    
    f=freq*2;
    % pad to remove influence of boundaries of padded series on the interior
    npad=3*round(1/f);
    nn=length(y);
    padded(npad+1:npad+nn)=y;
    padded([1:npad npad+nn+1:nn+2*npad])=0;
    [b,a]=butter(N,f,'low');
    yf0=filtfilt(b,a,padded);
    y_filtered=yf0(npad+1:nn+npad)';
    
elseif N<0 % Chebyshev method
    
    rip = 0.05;				% passband ripple in db
    [b,a] = cheby1(-N, rip, 2*freq);
    y_filtered = filtfilt(b, a, y);
    
    
end

% add back mean of unfiltered data
y_filtered=y_filtered+ym;
