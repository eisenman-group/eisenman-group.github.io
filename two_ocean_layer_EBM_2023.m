% This code numerically solves the two layer energy balance model (EBM)  
% described in Beer, Eisenman, Wagner and Fine (2023, hereafter BEWF23; see 
% reference below). This code was adapted from the code EBM_fast_WE15.m 
% which numerically solves the model described in Sec. 2b of Wagner & 
% Eisenman (2015; see reference below). Changes to the code include 
% removing an albedo change with surface temperture, adding a deeper ocean 
% layer to the model, and changing parameter values.
%
% This code allows some parameter values to be input. All other parameter 
% values are from BEWF23 (Table 1). The parameters to be input are: 
%
% F, which is a vector in years representing the radiative forcing (W m^-2),
% kv_w, vertical heat flux coefficient when T > 0 (W m^-2 K^-1),
% kv_i, vertical heat flux coefficient when T < 0 (W m^-2 K^-1),
% Ds, diffusion coefficient for heat transport in mixed layer (W m^-2 K^-1)
% Dd, diffusion coefficient for heat transport in deep layer (W m^-2 K^-1)
% A,  OLR when T = 0 (W m^-2)
%
% This code runs a simulation for the number of input values in vector F in 
% years with 5 years/timestep using a spatial resolution of 800 gridboxes,
% equally spaced in x=sin(lat), between the equator and pole.
%
% Emma Beer (ebeer@ucsd.edu), Adapted from EBM_fast_WE15.m Nov 2022.
%
% References: 
% E. Beer, I. Eisenman, T.J.W. Wagner and E.C. Fine (2023). A possible 
% hysteresis in the Arctic Ocean due to release of subsurface heat during 
% sea ice retreat. Journal of Physical Oceanography 53, 1323-1335.
%
% T.J.W. Wagner and I. Eisenman (2015). How climate model complexity 
% influences sea ice stability. Journal of Climate.

function [x_out,t_out,T_out,Td_out,Fb_out] = two_ocean_layer_EBM_2023(F_in,kv_w,kv_i,Ds,Dd,A)
%%Model parameters (BEWF23, Table 1) --------------------------------------
B  = 2.1;           % OLR temperature dependence (W m^-2 K^-1)
cw = 9.8*(50/75);   % ocean mixed layer heat capacity (W yr m^-2 K^-1), depth 50m
cwd = 9.8*(600/75); % ocean deep layer heat capacity (W yr m^-2 K^-1), depth 600m
S0 = 420;           % insolation at equator  (W m^-2)
S2 = 240;           % insolation spatial dependence (W m^-2)
a0 = 0.7;           % ice-free co-albedo at equator
a2 = 0.1;           % ice=free co-albedo spatial dependence
Tf = -2;
% -------------------------------------------------------------------------
n  = 800;           % grid resolution (number of points between equator and pole)
nt = 0.2;           % time resolution (time steps per year) - can be small since time integration is implicit
dur= length(F_in);
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;               %grid box width
x = (dx/2:dx:1-dx/2)';  %grid
%%Diffusion Operator (BEWF23, Appendix A) ---------------------------------
xb = (dx:dx:1.0-dx)';  
lambda=Ds/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);
%%Diffusion Operator for deep layer (BEWF23, Appendix A) ------------------
lambdad=Dd/dx^2*(1-xb.^2); L1d=[0; -lambdad]; L2d=[-lambdad; 0]; L3d=-L1d-L2d;
diffopd = - diag(L3d) - diag(L2d(1:n-1),1) - diag(L1d(2:n),-1);
% -------------------------------------------------------------------------
S = S0-S2*x.^2;         % insolation
T = -29*x.^2 -x + 24;   % IC for surface temperature
Td = -15*x.^2 -2*x + 20;% IC for deep layer temperature
aw = a0-a2*x.^2;        % co-albedo for open water
% -------------------------------------------------------------------------
allT = zeros(dur*nt,n);allTd = allT;
t = linspace(0,dur,dur*nt);
% integration over time using implicit difference and
% over x using central difference (through diffop and diffopd)
for i = 1:dur*nt
    F = F_in(floor((i-1)/nt)+1);
    kv = kv_w.*(T>Tf) + kv_i.*(T<=Tf);
    Fb = kv.*(Td-Tf) - kv.*(T-Tf).*(T>Tf);
    allFb(i,:) = Fb;
    % ---------------------------------------------------------------------
    % Rewrite equations (1) and (4) in BEWF23 in the format:
    % AA*T(n+1) = BB*T(n) + CC*Td(n+1) + DD
    % EE*Td(n+1) = FF*Td(n) + GG*T(n+1)
    I = eye(n);
    AA = I + dt/cw*(B*I - diffop + kv_w.*I.*(T>Tf));
    BB = I;
    CC = dt/cw*kv.*I;
    DD = dt/cw*(aw.*S - A + F);
    EE = I + dt/cwd*(-diffopd + kv.*I);
    FF = I;
    GG = dt/cwd*kv_w.*I.*(T>Tf);
    % ---------------------------------------------------------------------
    T = Tf + (AA - CC*inv(EE)*GG)\(BB*(T-Tf) + DD + CC*inv(EE)*FF*(Td-Tf));
    Td = Tf + inv(EE)*FF*(Td-Tf) + inv(EE)*GG*(T-Tf);
    allT(i,:)=T;
    allTd(i,:)=Td;
    % ---------------------------------------------------------------------
    if mod(i/nt,100)==0, disp(['year ' num2str(i/nt) ' complete']), end
end
% save results as output --------------------------------------------------
x_out=x';
T_out=allT;
Td_out=allTd;
Fb_out = allFb;
t_out = t;
