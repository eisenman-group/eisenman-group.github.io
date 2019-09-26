% EBM_fast_WE15.m:
% This code describes the standard annual-mean EBM described in Sec. 2b of
% the article referenced below, hereafter WE15. Here we use central
% difference spatial integration and Implicit Euler time stepping.
%
% The code EBM_simple_WE15.m, by contrast, uses a simpler formulation
% of the diffusion operator and time stepping with Matlab's ode45.
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
%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
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
nt = 5;       % time resolution (time steps per year) - can be small since time integration is implicit
dur= 30;
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;               %grid box width
x = (dx/2:dx:1-dx/2)';  %grid
%%Diffusion Operator (WE15, Appendix A) -----------------------------------
xb = (dx:dx:1.0-dx)';  
lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);

S = S0-S2*x.^2;        % insolation: WE15 eq. (3) with S_1 = 0
T = 10*ones(n,1);      % initial condition (constant temp. 10C everywhere)
aw = a0-a2*x.^2;       % co-albedo for open water

allT = zeros(dur*nt,n);
t = linspace(0,dur,dur*nt);
% integration over time using implicit difference and
% over x using central difference (through diffop)
for i = 1:dur*nt
    a = aw.*(T>0)+ai*(T<0); % WE15 eq. (4)
    C = a.*S-A+F;                    
    % Governing equation [cf. WE15 eq. (2)]:
    %    T(n+1) = T(n) + dt*(dT(n+1)/dt), with c_w*dT/dt=(C-B*T+diffop*T)
    % -> T(n+1) = T(n) + dt/cw*[C-B*T(n+1)+diff_op*T(n+1)]
    % -> T(n+1) = inv[1+dt/cw*(1+B-diff_op)]*(T(n)+dt/cw*C)
    T = (eye(n) + dt/cw*(B*eye(n)-diffop))\(T+dt/cw*C);
    allT(i,:)=T;
end
% plot results ------------------------------------------------------------
figure(1), clf
subplot(1,2,1)
plot(t,allT)
xlabel('t (years)'); 
ylabel('T (^oC, for different latitudes)')
xlim([0 t(end)])
subplot(1,2,2)
plot(x,T, 's-')
xlabel('x = sin \theta');
ylabel('T (^oC)')
box on
