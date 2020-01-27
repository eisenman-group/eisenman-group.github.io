% This code numerically solves the ocean-sea ice-climate model described
% in Beer, Eisenman, and Wagner (2020, hereafter BEW20; see reference 
% below) in the locked case with a vertical heat flux specified 
% by the input field Fb0. The case without locking is described in the code 
% sea_ice_EBM_deep_BEW20.m. This code was adapted from the code 
% sea_ice_EBM_WE15.m which numerically solves the model described in Wagner 
% & Eisenman (2015, hereafter WE15; see reference below). Changes to the 
% code include adding a deeper ocean layer to the model, changing parameter 
% values, and modifying the initial conditions.
%
% This code uses the default parameter values from BEW20 (Table S1). The
% two parameters to be input are:
%
% F, which is a single value representing the radiative forcing (W m^-2). 
% The default value is zero. 
%
% Fb0, which is the specified vertical heat flux field, with positive 
% defined as upwards, of size (n,nt). The default field used in BEW20 is 
% output as Fb1000_out from the code sea_ice_EBM_deep_BEW20.m when it is 
% run with no forcing (F=0).
%
% For computational efficiency, a diffusive 'ghost layer' is invoked, 
% as in WE15, and described in WE15 Appendix A. This allows us to solve the 
% system using an Implicit Euler timestep in the ghost layer (where
% diffusion occurs) and a Forward Euler timestep for the evolution of the 
% surface and deep layer enthalpy. The document WE15_NumericIntegration.pdf 
% provides further detailed documentation to go with the script 
% sea_ice_EBM_WE15.m, from which this code was adapted.
%
% By default this code runs a simulation for 500 years with 1000 
% timesteps/year using a spatial resolution of 100 gridboxes, equally 
% spaced in x=sin(lat), between the equator and pole.
%
% Emma Beer (ebeer@ucsd.edu), Adapted from sea_ice_EBM_WE15.m Jan 2020.
%
% References: 
% E. Beer, I. Eisenman, and T.J.W. Wagner (2020). Polar amplification due 
% to enhanced heat flux across the halocline. Geophys Res Lett.
%
% T.J.W. Wagner and I. Eisenman (2015). How climate model complexity 
% influences sea ice stability. J Climate.
%
%--------------------------------------------------------------------------
function [tfin_out, Efin_out, Edfin_out, Tfin_out, Tdfin_out, Fbfin_out] = sea_ice_EBM_deep_locked_BEW20(F,Fb0)
%%Model parameters (BEW20 Table S1) -------------------------
D  = 0.5;     %diffusivity for heat transport (W m^-2 K^-1)
A  = 189;     %OLR when T = T_m (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8*(50/75);     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420;     %insolation at equator  (W m^-2)
S1 = 354.9;     %insolation seasonal dependence (W m^-2)
S2 = 240;     %insolation spatial dependence (W m^-2)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice-free co-albedo spatial dependence
ai = 0.4;     %co-albedo where there is sea ice
k  = 2;       %sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)
cg = 0.01*cw; %ghost layer heat capacity(W yr m^-2 K^-1)
tau = 1e-5;   %ghost layer coupling timescale (yr)
Dd  = 0.08;     %diffusivity for deep heat transport (W m^-2 K^-1)
cwd = 9.8*(600/75);  %ocean mixed layer heat capacity for layer depth 600 (W yr m^-2 K^-1)
%kv = 2; %vertical heat flux coefficient (W m^-2 K^-1)
cgd = 0.01*cwd; %deep ghost layer heat capacity(W yr m^-2 K^-1)
Tm = -2; %melting point for Arctic conditions (K)

%%The default run in BEW20 Fig. 2 uses these time-stepping parameters: ----
% n  = 400;   % # of evenly spaced latitudinal gridboxes (equator to pole)
% nt = 1e3;   % # of timesteps per year (limited by numerical stability)
% dur= 2000;   % # of years for the whole run
%%For a quicker computation, here we use these parameters: ----------------
n  = 100;   % # of evenly spaced latitudinal gridboxes (equator to pole)
nt = 1e3;   % # of timesteps per year (limited by numerical stability)
dur= 500;    % # of years for the whole run
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;               %grid box width
x = (dx/2:dx:1-dx/2)';  %grid
%%Diffusion Operator (WE15 Appendix A) ------------------------------------
xb = (dx:dx:1.0-dx)';  
lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);

%%Deep layer Diffusion Operator -------------------------------------------
lambda_d=Dd/dx^2*(1-xb.^2); L1d=[0; -lambda_d]; L2d=[-lambda_d; 0]; L3d=-L1d-L2d;
diffop_d = - diag(L3d) - diag(L2d(1:n-1),1) - diag(L1d(2:n),-1);

%%Definitions for implicit scheme for Tg
cg_tau = cg/tau;
dt_tau = dt/tau;
dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;

%%Definitions for implicit scheme for Tgd
cg_tau_d = cgd/tau;
kappa_d = (1+dt_tau)*eye(n)-dt*diffop_d/cgd;

%%Seasonal forcing [BEW20 Eq. (10), WE15 Eq. (3)]
ty = dt/2:dt:1-dt/2;
S=repmat(S0-S2*x.^2,[1,nt])-repmat(S1*cos(2*pi*ty),[n,1]).*repmat(x,[1,nt]);
%%Further definitions
M = B+cg_tau;
Md = cg_tau_d;
aw= a0-a2*x.^2;   %open water albedo
kLf = k*Lf;
%%Initial conditions ------------------------------------------------------
T = 25.5+19.9.*x-65.9.*x.^2 + F/B;
Td = 25.9-5.5.*x-22.3.*x.^2 + F/B;
Tg = T; E = cw*(T-Tm);Tgd = Td; Ed = cwd*(Td-Tm);
%%Set up output arrays, saving 100 timesteps/year
E100 = zeros(n,dur*100); T100 = E100;
Ed100 = E100; Td100 = Ed100; Fb100 = E100;
%%Integration (see WE15_NumericIntegration.pdf)----------------------------
% Loop over Years ---------------------------------------------------------
for years = 1:dur
    % Loop within One Year-------------------------------------------------
    for i = 1:nt
        if mod(i,nt/100)==0 %store 100 timesteps per year
            E100(:,(years-1)*100+i/(nt/100)) = E;
            T100(:,(years-1)*100+i/(nt/100)) = T;
            Ed100(:,(years-1)*100+i/(nt/100)) = Ed;
            Td100(:,(years-1)*100+i/(nt/100)) = Td;
            Fb100(:,(years-1)*100+i/(nt/100)) = Fb0(:,i);
        end
        % forcing
        alpha = aw.*(E>0) + ai*(E<0);       %BEW20 Eq. (9), WE15 Eq. (4)
        C = alpha.*S(:,i) + cg_tau*Tg-A+F+B*Tm; 
        Cd = cg_tau_d*Tgd; 
        % surface temperature
        T0 =  (C-kLf.*Tm./E)./(M-kLf./E);                 %WE15 Eq. (A3)
        T = (Tm+E/cw).*(E>=0)+Tm.*(E<0).*(T0>Tm)+T0.*(E<0).*(T0<Tm); %BEW20 Eq. (7), WE15 Eq. (9)
        Td = (Tm+Ed/cwd);
        % Forward Euler for E and Ed. Here, the vertical heat flux field 
        % (Fb0) is used instead of a flux proportional to the temperature 
        % of the two layers.
        E = E+dt*(C-M*T+Fb0(:,i)); 
        Ed = Ed+dt*(Cd-Md*Td-Fb0(:,i)); 
        % Implicit Euler for Tg and Tgd
        Tg = (kappa-diag(dc./(M-kLf./E).*(T0<Tm).*(E<0)))\ ...
            (Tg + (dt_tau*((Tm+E/cw).*(E>=0)+(ai*S(:,i) ...
            -A+F+B*Tm-kLf.*Tm./E)./(M-kLf./E).*(T0<Tm).*(E<0)+Tm.*(E<0).*(T0>Tm))));        %WE15 Eq. (A1)
        Tgd = kappa_d\(Tgd + dt_tau*(Tm+Ed/cwd));
    end
    if mod(years,10)==0, disp(['year ' num2str(years) ' complete']), end
end
% -------------------------------------------------------------------------
%%output only final year 
tfin = linspace(0,1,100);  
Efin = E100(:,end-99:end); 
Tfin = T100(:,end-99:end);
Edfin = Ed100(:,end-99:end); 
Tdfin = Td100(:,end-99:end);
Fbfin = Fb100(:,end-99:end);
if nargout>0 %save output
    tfin_out=tfin;
    Efin_out=Efin;
    Edfin_out=Edfin;
    Tfin_out=Tfin;
    Tdfin_out=Tdfin;
    Fbfin_out=Fbfin;
end
% -------------------------------------------------------------------------
