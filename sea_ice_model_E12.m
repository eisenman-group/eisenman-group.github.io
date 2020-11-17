function [final_year, full_time_series]=sea_ice_model_E12(varargin)
%   [final_year, full_time_series]=sea_ice_model_E12(arguments)
%
% Toy model of sea ice in an ocean mixed layer below a thermodynamic
% atmosphere. Nondimensional parameters and model state are described in
% the reference below.
%
%   Output: final_year or full_time_series = [t E Tsrf]
%
% Example: Plot steady-state seasonal cycle with seasonally ice-free conditions
%   Y=sea_ice_model_E12('Lm=0.95','FB=0');
%   subplot(2,1,1), plot(Y(:,1),-Y(:,2)*5.9), ylabel('ice thickness (m)')
%   subplot(2,1,2), plot(Y(:,1),Y(:,3)*8.8), ylabel('sfc temp (^oC)')
%
% See also the dimensional version of this model at
%   http://eisenman.ucsd.edu/code/sea_ice_model_E12_dimensional.m
%
% Reference:
%   Ian Eisenman, 2012. Factors controlling the bifurcation structure of
%   sea ice retreat. J Geophys Res 117, D01111. 
%
% Ian Eisenman (eisenman@ucsd.edu), 2012

global Sa Lm La phi B zeta Da ha FB

% === Model parameters ===
Sa=1.5;
Lm=1.25;
La=0.73;
phi=0.15;
B=0.45;
zeta=0.12;
Da=0.43;
ha=0.08;
FB=0;

% === Integration parameters ===
E0=-0.5; % initial value of E
max_dur=50; % max duration of integration (yrs)
spin_up=1; % whether to stop before max_dur if model is spun up
spun_up=3e-4; % tolerance for model to be considered spun up
silent=0; % if 1, don't display error or final value
odesolver='ode45'; % which ode solver to use
RelTol=1e-5; % for ODE solver
AbsTol=1e-6; % for ODE solver
t0=0; % time offset for starting point

% === EDIT PARAMETERS BELOW HERE ===
% parameter value changes as input to function (batch mode)
if nargin>0, for j=1:length(varargin), eval([varargin{j} ';']), end, end
% parameter value changes for interactive run
if nargout==0
    % = Enter changes to default parameter values here for interactive run =
    %Lm=1.1;
end

% ======

% === Integration ===
% run integration one yr at a time until spun_up condition or max duration
intyr=1; E=[]; t=[];
options=odeset('RelTol',RelTol,'AbsTol',AbsTol);
while intyr<=max_dur
    if strcmp(odesolver,'ode45'), [tt,EE]=ode45(@dEdt,intyr+t0+[-1 0],E0,options); end
    if strcmp(odesolver,'ode23t'), [tt,EE]=ode23t(@dEdt,intyr+t0+[-1 0],E0,options); end
    if strcmp(odesolver,'ode15s'), [tt,EE]=ode15s(@dEdt,intyr+t0+[-1 0],E0,options); end
    if strcmp(odesolver,'ode23'), [tt,EE]=ode23(@dEdt,intyr+t0+[-1 0],E0,options); end
    E=[E; EE]; t=[t; tt];
    intyr2=intyr;
    E0=EE(end,:); intyr=intyr+1;
    if ( abs(EE(end,:)-EE(1,:)) < spun_up ) && spin_up, intyr=max_dur+1; end
end
if silent==0 % display final value/error in command window
    if intyr2<max_dur
        disp(['Integration ended at year ' num2str(intyr2) ' with final value E0=' mat2str(EE(end),2)])
    else
        disp(['Integration ended at year ' num2str(intyr2) ' with final value E0=' mat2str(EE(end),2)])
        disp(['ERROR: Integral did not converge during ' num2str(intyr2) ' years.'])
    end
end

% === plotting or output ===
% find sfc temp from t,E output using model function
[F TT]=dEdt(tt,EE);
tt=tt-tt(1);
Yf=[tt EE TT];
[F T]=dEdt(t,E);
Y=[t E T];
if nargout>0 % output
    final_year=Yf;
    full_time_series=Y;
else % plotting
    save tmp.mat Yf Y
    % plot evolution of model state: full time series
    if sum(2==get(0,'children')), close(2), end
    figure(2), clf, set(gcf,'position',[19 407 560 402])
    subplot(2,2,1)
    plot(t,E), title('E')
    axis tight, grid on
    subplot(2,2,2)
    plot(t,T), title('Tsrf')
    axis tight, grid on
    % plot evolution of model state: final year
    subplot(2,2,3)
    plot(tt,EE), title('E')
    axis tight, grid on
    set(gca,'xtick',(15:30:360)/360,'xticklabel',[])
    subplot(2,2,4)
    TT(EE>=0)=NaN;
    plot(tt,TT), title('Tsrf')
    axis tight, grid on, yl=get(gca,'ylim'); if yl(2)<0, yl(2)=0; end, set(gca,'ylim',yl)
    set(gca,'xtick',(15:30:360)/360,'xticklabel',[])
end


% ===================


% === model equations ===
function [F T]=dEdt(t,E)
global Sa Lm La phi B zeta Da ha FB

A=(1+Da*tanh(E/ha)).*(1-Sa*cos(2*pi*t))+(-Lm-La*cos(2*pi*(t-phi)));
T= (E>=0).*E + (E<0).*(A<0).*( A./B.*(1-zeta./E).^(-1) );
F=A-B.*T+FB;
