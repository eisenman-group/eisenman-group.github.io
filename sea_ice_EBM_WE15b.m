% This script numerically solves the diffusive energy balance model (EBM)
% with seasonal variations, sea ice, and stochastic weather noise described
% in eqn (1) of Wagner and Eisenman (2015b).
%
% The original version of this model was introduced and described in Wagner
% and Eisenman (2015a, hereafter WE15). This version of the model has the
% addition of stochastic weather noise.
%
% Till Wagner (tjwagner@ucsd.edu) & Ian Eisenman (eisenman@ucsd.edu), 
% created Nov 2015, minor bug fix Jan 2022 [in Eq. (A1), S(:,i) -> S(:,i+1)].
%
% References:
% T.J.W. Wagner and I. Eisenman (2015a). How climate model complexity
%   influences sea ice stability. J Climate 28, 3998-4014.
% T.J.W. Wagner and I. Eisenman (2015b). False alarms: How early warning
%   signals falsely predict abrupt sea ice loss. Geophys Res Lett 42, 10333-10341.
 
n   = 200;   %grid resolution
dur = 350;   %duration of simulation
sig = 0.5;   %noise amplitude
Fdef= 0;     %initial Forcing level 
             %(F = 0 corresponds roughly to pre-industrial levels)
spinup = 50; %start ramping after 'spinup' years
%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
D  = 0.6;     %diffusivity for heat transport (W m^-2 K^-1)
S1 = 338;     %insolation seasonal dependence (W m^-2)
A  = 193;     %OLR when T = T_m (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420;     %insolation at equator  (W m^-2)
S2 = 240;     %insolation spatial dependence (W m^-2)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice=free co-albedo spatial dependence
ai = 0.4;     %co-albedo where there is sea ice
Fb = 4;       %heat flux from ocean below (W m^-2)
k  = 2;       %sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)
cg = 1e-3;    %ghost layer heat capacity(W yr m^-2 K^-1)
tau = 3e-6;   %ghost layer coupling timescale (yr)
%%time stepping
nt = 1e3;
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;  x = (dx/2:dx:1-dx/2)';  %native grid
%%Diffusion Operator (WE15, Appendix A) -----------------------------------
xb = (dx:dx:1.0-dx)';  lambda=D/dx/dx*(1-xb.*xb);
a=[0; -lambda]; c=[-lambda; 0]; b=-a-c;
diffop = - diag(b) - diag(c(1:n-1),1) - diag(a(2:n),-1);
%%Definitions for implicit scheme on Tg
cg_tau = cg/tau; dt_tau = dt/tau; dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;
%%Seasonal forcing (WE15 eq.3)
ty = dt/2:dt:1-dt/2;
S=repmat(S0-S2*x.^2,[1,nt])-repmat(S1*cos(2*pi*ty),[n,1]).*repmat(x,[1,nt]);
S=[S S(:,1)];
%%Further definitions
M = B+cg_tau; 
aw= a0-a2*x.^2; 
kLf = k*Lf;
%%Set up output arrays, saving 100 timesteps/year
E100 = zeros(n,dur*100); 
T100 = zeros(n,dur*100); 
%%ramping rate
dF = 1/(20*nt);
%%Initial conditions ------------------------------------------------------
T = 10*ones(n,1);
Tg = T; E = cw*T;
%%noise timeseries
sig_noise = sig/sqrt(dt);
noise = sig_noise*randn(1,dur*nt);
lp = 1/52;  %1 week 'decorrelation' time
alpha = exp(-dt/lp); nalpha = sqrt(1-alpha^2);
N_red = noise*0;
N_red(1) = noise(1);
for i=2:length(noise)
    N_red(i)=alpha*N_red(i-1)+nalpha*noise(i);
end
%%run the model -----------------------------------------------------------
p = 0; m = 0; N = 0; F = Fdef;
%%Integration (see WE15_NumericIntegration.pdf)----------------------------
%%Loop over Years ---------------------------------------------------------
tic
for years = 1:dur
    %Loop within One Year-------------------------------------------------
    for i = 1:nt
        m = m+1;
        if years > spinup  %start ramping and noise after spin up
            F = F+dF;
            N = N_red(m);
        end
        %store 100 timesteps per year
        if (p+1)*nt/100 == m
            p = p+1;
            E100(:,p) = E;
            T100(:,p) = T;
        end
        % forcing
        alpha = aw.*(E>0) + ai*(E<0);        %WE15, eq.4
        C =alpha.*S(:,i)+cg_tau*Tg-A+F+N;
        % surface temperature
        T0 =  C./(M-kLf./E);                 %WE15, eq.A3
        T = E/cw.*(E>=0)+T0.*(E<0).*(T0<0);  %WE15, eq.9
        % Forward Euler on E
        E = E+dt*(C-M*T+Fb);                 %WE15, eq.A2
        % Implicit Euler on Tg
        Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\ ...
            (Tg + (dt_tau*(E/cw.*(E>=0)+(ai*S(:,i+1) ...
            -A+F+N)./(M-kLf./E).*(T0<0).*(E<0))));        %WE15, eq.A1
    end
    yrs = sprintf('year %d complete',years); disp(yrs)
end
toc
%%--------------------------------------------------------------------------
%compute ice edge
xi = ones(1,p);
for i = 1:p
    if any(E100(:,i)<0)==1
    xi(i) = x(find(E100(:,i)<0,1,'first'));
    end
end
%--------------------------------------------------------------------------
%compute yearly summer and winter ice area and temperature at pole
SIA = zeros(1,dur-spinup);
Tpole = SIA;
        sept= (76+spinup*100):100:p; %summer dates excluding spinup
        mar = (26+spinup*100):100:p; %winter dates excluding spinup 
        SIA_sept= 255*(1-xi(sept));  %summer ice area in million km^2
        SIA_mar = 255*(1-xi(mar));   %winter ice area in million km^2
        Tpole_sept= T100(end,sept);  %summer temperature at pole
        Tpole_mar = T100(end,mar);   %winter temperature at pole
 
Fv = linspace(Fdef,F,dur-spinup);   %yearly forcing without spinup
%--------------------------------------------------------------------------
%plot summer/winter ice areas and temperatures at the pole
figure(1); clf
subplot(2,1,1)
plot(Fv, SIA_sept,'r',Fv,SIA_mar,'b')
xlabel('F (W m^{-2)}')
ylabel('Sea Ice Area (10^6 km^2)')
legend('summer','winter')
 
subplot(2,1,2)
plot(Fv, Tpole_sept,'r',Fv,Tpole_mar,'b')
xlabel('F (W m^{-2)}')
ylabel('Temperature at Pole (^oC)')
