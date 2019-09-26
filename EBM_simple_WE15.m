% EBM_simple_WE15.m:
% This code describes the standard annual-mean EBM described in Sec. 2b of
% the article referenced below, hereafter WE15. Here we use central
% difference spatial integration and time stepping with MATLAB's ode45.
%
% The code EBM_fast_WE15.m, by contrast, uses a faster but more complicated
% formulation of the diffusion operator and Implicit Euler time stepping.
%
% Parameters are as described in WE15, Table 1. Note that we do not include 
% ocean heat flux convergence or a seasonal cylce in the forcing 
% (equivalent to S_1 = F_b = 0 in WE15). This code uses an albedo for ice
% when T<0 (WE15 instead uses the condition E<0, since it includes a
% representation of ice thickness). In this code, we define T = Ts - Tm,
% where Ts is the surface temperature and Tm the melting point (WE15, by
% contrast, defines T = Ts).
%
% Till Wagner & Ian Eisenman, Mar 2015
% tjwagner@ucsd.edu or eisenman@ucsd.edu
%
% Reference: "How Model Complexity Influences Sea Ice Stability", 
% T.J.W. Wagner & I. Eisenman, J Clim (2015).
%
function [x_out,T_out] = EBM_simple_WE15
%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
global  D B n cw aw ai S A F x dx 
D  = 0.6;     % diffusivity for heat transport (W m^-2 K^-1)
A  = 193;     % OLR when T = 0 (W m^-2)
B  = 2.1;     % OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     % ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420;     % insolation at equator  (W m^-2)
S2 = 240;     % insolation spatial dependence (W m^-2)
a0 = 0.7;     % ice-free co-albedo at equator
a2 = 0.1;     % ice=free co-albedo spatial dependence
ai = 0.4;     % co-albedo where there is sea ice
F  = 0;       % radiative forcing (W m^-2)
% -------------------------------------------------------------------------
n  = 50;      % grid resolution (number of points between equator and pole)
x  = linspace(0,1,n);
dx = 1/(n-1);
S  = S0-S2*x.^2;      % insolation [WE15 eq. (3) with S_1 = 0]
aw = a0-a2*x.^2;      % open water co-albedo
T0 = 10*ones(1,n);    % initial condition (constant temp. 10C everywhere)
tspan = [0 30];       % timespan in years [t_start t_end]
% time integration---------------------------------------------------------
[t,T] = ode45(@rhsWE15,tspan,T0);
if nargout==0 % plot results-----------------------------------------------
    figure(1), clf
    subplot(1,2,1)
    plot(t,T)
    xlabel('t (years)'); 
    ylabel('T (^oC, for different latitudes)')
    xlim([0 t(end)])
    box on
    subplot(1,2,2)
    hold all
    plot(x,T(end,:), 's-')
    xlabel('x = sin \theta');
    ylabel('T (^oC)')
    box on
else % save results--------------------------------------------------------
    x_out=x;
    T_out=T;
end
% -------------------------------------------------------------------------
% ODE with spatial finite differencing for ode45 --------------------------
function dTdt = rhsWE15(~,T)
global  D B n cw aw ai S A F x dx
Tdot = zeros(n,1);
% forcing
alpha = aw.*(T>0)'+ai*(T<0)';   % WE15 eq. (4)
C = alpha.*S-A+F;
% solve c_w dT/dt = D (1-x^2) d^2 T/dx^2 - 2 x D dT/dx + C - B T [cf. WE15 eq. (2)]
% use central difference
Tdot(2:n-1) = (D/dx^2)*(1-x(2:n-1).^2)'.*(T(3:n)-2*T(2:n-1)+T(1:n-2)) ...
              -(D*x(2:n-1)/dx)'.*(T(3:n)-T(1:n-2));
% boundary condition: let dT/dt = 0 at x=0
Tdot(1) = D*(2*(T(2)-T(1))/dx^2);   % WE15 eq. (7)
% use asymmetric (backward) spatial difference for the pole
Tdot(n) = - D*(2*x(n)*(T(n)-T(n-1))/dx);    %O(dt) accuracy at x=1
% temperature tendency
dTdt = (Tdot+C' - B*T)/cw;
