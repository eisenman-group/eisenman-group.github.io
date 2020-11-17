% Dimensional version of the sea ice toy model at
%   http://eisenman.ucsd.edu/code/sea_ice_model_E12.m
%
% Reference:
%   Ian Eisenman, 2012. Factors controlling the bifurcation structure of
%   sea ice retreat. J Geophys Res 117, D01111. 
%
% Ian Eisenman (eisenman@ucsd.edu), 2020

a_bar = 0.56;   % coalbedo averaged between ice and ocean
delta_a = 0.48; % difference in coalbedo between ice and ocean
h_alpha = 0.5;	% smoothness of albedo transition (m)
B = 2.83;       % dependence of net surface flux on surface temperature (W/m^2/K)
Fb = 0;         % upward heat flux into bottom
Sm = 100;       % downward shortwave annual-mean (W/m^2)
Sa = 150;       % downward shortwave seasonal amplitude (W/m^2)
Lm = 70;        % reference longwave annual-mean (W/m^2)
La = 41;        % reference longwave seasonal amplitude (W/m^2)
P = 1;          % forcing period (yrs)
phi = 0.15;     % time between summer solstice and peak of longwave forcing (yrs)
coHo = 6.3;     % heat capacity of ocean mixed layer (W/m^2 yr/K) (equivalent to 2e8 J/m^2/K)
Li = 9.5;       % sea ice latent heat of fusion (W/m^3 yr) (equivalent to 3e8 J/m^3)
zeta = 0.7;     % sea ice thermodynamic scale thickness zeta=ki/B (m)
E0=-26.16;      % initial condition

A = @(t,E)( (a_bar + delta_a/2 * tanh(E/(Li*h_alpha)) ).*( Sm - Sa.*cos(2*pi*t/P) ) ...
    - (Lm +La.*cos(2*pi*(t/P-phi/P)) ) ); % E12 eq. (4)
T = @(t,E)( E/coHo.*(E>=0) + A(t,E)/B.*(1-zeta.*Li./E).^(-1).*(E<0).*(A(t,E)<0) ); % E12 eq. (6) with hi=-E/Li from E12 eq. (1)
dEdt = @(t,E)( A(t,E) - B*T(t,E) + Fb ); % E12 eqs. (2)-(3)

[t, E] = ode45(dEdt,[0 1],E0,odeset('RelTol',1e-5,'AbsTol',1e-6));
hi=0-E/Li.*(E<0); % E12 eq. (1)

subplot(1,2,1)
plot(t,hi), ylabel('h_i (m)'), xlabel('t (yrs)')
subplot(1,2,2)
plot(t,T(t,E)), ylabel('T (^oC)'), xlabel('t (yrs)')
