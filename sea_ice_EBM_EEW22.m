% Reference: "Spurious Climate Impacts in Coupled Sea Ice Loss Simulations".
% M.R. England, I. Eisenman, and T.J.W. Wagner. J Clim (2022).
%
% Adapted from the model code described in "How Model Complexity Influences
% Sea Ice Stability", T.J.W. Wagner & I. Eisenman, J Clim (2015).
%
% The default configuration here runs a simulation for 200 years at 1000
% timesteps/year and a spatial resolution of 400 gridboxes, equally spaced
% between equator and pole.
%
% The sea_ice_EBM_EEW22 function takes in the following inputs:
% (experiment,F,Fb,dur,nt,n,dummyX)
%
% The experiment input describes to the type of experiment
% 0 = default WE15 simulations
% 1 = prescribed albedo simulations
% 2 = albedo modification method
% 3 = ghostflux method
% 4 = nudging method
%
% F = climate forcing (W m^-2)
% Fb = heat flux from deep ocean (W m^-2)
% dur = number of years to run simulation for
% nt = timesteps per year
% n = number of gridboxes equally spaced between equator and pole
%
% The dummy1 and dummy2 input variables depends on the experiment
% if experiment=0 -> dummy1 and dummy2 are not used
% if experiment=1 -> dummy1 is the enthalpy field E from the simulation you
%   want to specify the albedo from, dimensions [n, nt]
% if experiment=2 -> dummy1 is the value for the new ice coalbedo (1 -
%   albedo), which has a default value of 0.4
% if experiment=3 -> dummy1 is the time varying ghostflux, dimensions [nt]
%                 -> dummy2 is enthalpy field E from the simulation you
%                 want to simulate the sea ice from, dimensions [n, nt]
% if experiment=4 -> dummy1 is the enthalpy field E from the simulation you
%   want to nudge the sea ice to, dimensions [n, nt]
%                 -> dummy2 is the nudging timescale in terms of timesteps
%
% Output
% E: Enthalpy
% T: Temperature
% f: fluxes
% h: ice thickness
% (fin): final year at coarse temporal resolution
% (finnt): final year at high temporal resolution
%
% Mark England, 28th October 2022
% markengland20@gmail.com, markrossengland.com, @markrossengland (twitter)
%
%--------------------------------------------------------------------------

function [tfin, Efin, Efinnt, Tfin,Tfinnt,ffin,ffinnt,hfin,hfinnt] = sea_ice_EBM_EEW22(experiment,F,Fb,dur,nt,n,dummy1,dummy2)

%%Model parameters (WE15, Table 1 and Section 2d) -------------------------

D  = 0.6;     %diffusivity for heat transport (W m^-2 K^-1)

A  = 193.0-0.35;     %OLR when T = T_m (W m^-2)

B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)

cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)

S0 = 420;     %insolation at equator  (W m^-2)

S1 = 338;     %insolation seasonal dependence (W m^-2)

S2 = 240;     %insolation spatial dependence (W m^-2)

a0 = 0.7;     %ice-free co-albedo at equator

a2 = 0.1;     %ice-free co-albedo spatial dependence

ai = 0.4;     %co-albedo where there is sea ice

% Fb = 4;       %heat flux from ocean below (W m^-2)

k  = 2;       %sea ice thermal conductivity (W m^-2 K^-1)

Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)

% Tm = 0;     %melting temp., not included in the eqns below

% F  = 0;       %radiative forcing (W m^-2)

cg = 0.01*cw; %ghost layer heat capacity(W yr m^-2 K^-1)

tau = 1e-5;   %ghost layer coupling timescale (yr)

%%The default run in WE15, Fig 2 uses the time-stepping parameters: -------

% n=400; % # of evenly spaced latitudinal gridboxes (equator to pole)

% nt=1e3; % # of timesteps per year (approx lower limit of stability)

% dur=200; % # of years for the whole run

%%For a quicker computation, use the parameters: --------------------------

% n  = 100;
% 
% nt = 1e3;
% 
% dur= 30;

dt = 1/nt;

%%Spatial Grid ------------------------------------------------------------

dx = 1/n;     %grid box width

x = (dx/2:dx:1-dx/2)';  %native grid

E_target_min = -5; % ghost flux will be still applied in areas of remaining thin sea ice (W m^-2 K^-1)

%%Diffusion Operator (WE15, Appendix A) -----------------------------------

xb = (dx:dx:1.0-dx)';  

lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;

diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);

%%Definitions for implicit scheme on Tg

cg_tau = cg/tau;

dt_tau = dt/tau;

dc = dt_tau*cg_tau;

kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;

%%Seasonal forcing (WE15 eq.3)

ty = dt/2:dt:1-dt/2;

S=repmat(S0-S2*x.^2,[1,nt])-repmat(S1*cos(2*pi*ty),[n,1]).*repmat(x,[1,nt]);

%%Further definitions

M = B+cg_tau;

aw= a0-a2*x.^2;   %open water albedo

kLf = k*Lf;

%%Initial conditions ------------------------------------------------------

T = 7.5+20*(1-2*x.^2);

Tg = T; E = cw*T; C= zeros(size(T));

%%Set up output arrays, saving 100 timesteps/year

vec100 = 1:nt/100:nt*dur;

E100 = zeros(n,dur*100); T100 = E100; f100 = E100;

vecnt = 1:1:nt*dur;

Ent = zeros(n,dur*nt); Tnt = Ent; fnt = Ent; 

p = 0; m = 0; q = 0;

%%Integration (see WE15_NumericIntegration.pdf)----------------------------

% Loop over Years ---------------------------------------------------------

for years = 1:dur

    % Loop within One Year-------------------------------------------------

    for i = 1:nt

        %store 100 timesteps per year

        m = m+1;

        if (p+1)*nt/100 == m

            p = p+1;

            E100(:,p) = E;

            T100(:,p) = T;
            
            f100(:,p) = C-M*T+Fb;

        end
      
        
        if (q+1) == m

            q = q+1;

            Ent(:,q) = E;

            Tnt(:,q) = T;
            
            fnt(:,q) = C-M*T+Fb;

        end

        if experiment==2
            ai = dummy1;
        end

        % forcing
        if experiment==1
            alpha = aw.*(dummy1(:,i)>0) + ai*(dummy1(:,i)<0);    % WE15, eq.4
        else
            alpha = aw.*(E>0) + ai*(E<0);    % WE15, eq.4
        end

        if experiment==3
            C =alpha.*S(:,i)+cg_tau*Tg-A+F+dummy1(:,i).*(E<0).*(dummy2(:,i)>E_target_min);  
        else
            C =alpha.*S(:,i)+cg_tau*Tg-A+F;
        end


        % surface temperature

        T0 =  C./(M-kLf./E);                 %WE15, eq.A3

        T = E/cw.*(E>=0)+T0.*(E<0).*(T0<0);  %WE15, eq.9

        % Forward Euler on E

        if experiment==4
            E = E+dt*(C-M*T+Fb)+(1/dummy2)*(1-E).*(E<0&dummy1(:,i)>0);                 %WE15, eq.A2
        else
            E = E+dt*(C-M*T+Fb);                 %WE15, eq.A2
        end

        

        % Implicit Euler on Tg

        if i==nt
            ii=0;
        else
            ii=i;
        end

        if experiment==3
            Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\ ...
                (Tg + (dt_tau*(E/cw.*(E>=0)+(ai*S(:,ii+1) ...
                -A+F+dummy1(:,i).*(E<0).*(dummy2(:,i)>E_target_min))./(M-kLf./E).*(T0<0).*(E<0))));        %WE15, eq.A1
        else
            Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\ ...
                (Tg + (dt_tau*(E/cw.*(E>=0)+(alpha.*S(:,ii+1) ...
                -A+F)./(M-kLf./E).*(T0<0).*(E<0))));        %WE15, eq.A1
        end

    end

%     yrs = sprintf('year %d complete',years); disp(yrs)

end
disp("done")
% -------------------------------------------------------------------------

%%output only converged, final year 

tfin = linspace(0,1,100);  
tfinnt = linspace(0,1,nt); 

Efin = E100(:,end-99:end); 
Efinnt = Ent(:,(end-nt+1):end); 

Tfin = T100(:,end-99:end);
Tfinnt = Tnt(:,(end-nt+1):end);

ffin = f100(:,end-99:end);
ffinnt = fnt(:,(end-nt+1):end);

hfin = -(Efin/Lf.*(Efin<0));
hfinnt = -(Efinnt/Lf.*(Efinnt<0));

 
