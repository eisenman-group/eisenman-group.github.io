%%
% This script showcases how to use the function sea_ice_EBM_EEW22 to
% replicate the published model results.
% Reference: "Spurious Climate Impacts in Coupled Sea Ice Loss Simulations". 
% M.R. England, I. Eisenman, and T.J.W. Wagner. J Clim (2022).

%%
clear all
clc
close all
% 
set(0,'defaultfigurecolor',[1 1 1])

% The default configuration here runs a simulation for 200 years at 1000
% timesteps/year and a spatial resolution of 400 gridboxes, equally spaced
% between equator and pole.
%
% The results of the study function takes in the following inputs:
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
% acronym guide
% df = default WE15 run
% al = modified albedo
% ng = nudging
% ht = ghost flux heating
%
% Mark England, 28th October 2022
%
% markengland20@gmail.com, markrossengland.com, @markrossengland (twitter)

%% Parameters

Fb=4;
dur = 200;
nt = 1e3;
n=400;
w=2*pi; % annual frequency
B=2.1; % sensitivity of OLR to sfc temp
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
D  = 0.6;     %diffusivity for heat transport (W m^-2 K^-1)

%% Default WE15 simulations

% Default simulation with 0 W/m2 forcing
[tfin, df_Efin_F0,df_Efinnt_F0,df_Tfin_F0,df_Tfinnt_F0,df_ffin_F0,df_ffinnt_F0,df_hfin_F0,df_hfinnt_F0] = sea_ice_EBM_EEW22(0,0.0,Fb,dur,nt,n,0,0);

% Default simulation with 3.1 W/m2 forcing
[tfin, df_Efin_F3,df_Efinnt_F3,df_Tfin_F3,df_Tfinnt_F3,df_ffin_F3,df_ffinnt_F3,df_hfin_F3,df_hfinnt_F3] = sea_ice_EBM_EEW22(0,3.1,Fb,dur,nt,n,0,0);

%% Specified albedo simulations

% Default WE15 with ice albedo specified from future
[tfin, df_Efin_F0_albedo3,df_Efinnt_F0_albedo3,df_Tfin_F0_albedo3,df_Tfinnt_F0_albedo3,df_ffin_F0_albedo3,df_ffinnt_F0_albedo3,df_hfin_F0_albedo3,df_hfinnt_F0_albedo3] = sea_ice_EBM_EEW22(1,0,Fb,dur,nt,n,df_Efinnt_F3,0);
% Default WE15 with ice albedo specified from control and forcing from
% future400
[tfin, df_Efin_F3_albedo0,df_Efinnt_F3_albedo0,df_Tfin_F3_albedo0,df_Tfinnt_F3_albedo0,df_ffin_F3_albedo0,df_ffinnt_F3_albedo0,df_hfin_F3_albedo0,df_hfinnt_F3_albedo0] = sea_ice_EBM_EEW22(1,3.1,Fb,dur,nt,n,df_Efinnt_F0,0);

%% Albedo modification method SEA ICE LOSS
% matches the summer minimum;
[tfin, al_Efin_F0_a45,al_Efinnt_F0_a45,al_Tfin_F0_a45,al_Tfinnt_F0_a45,al_ffin_F0_a45,al_ffinnt_F0_a45,al_hfin_F0_a45,al_hfinnt_F0_a45] = sea_ice_EBM_EEW22(2,0,Fb,dur,nt,n,0.45,0);
% matches the annual mean;
[tfin, al_Efin_F0_a48,al_Efinnt_F0_a48,al_Tfin_F0_a48,al_Tfinnt_F0_a48,al_ffin_F0_a48,al_ffinnt_F0_a48,al_hfin_F0_a48,al_hfinnt_F0_a48] = sea_ice_EBM_EEW22(2,0,Fb,dur,nt,n,0.48,0);

%% Ghostflux method SEA ICE LOSS
forcing_arctic = 35; % W/m2

A = 30;
B = 0.02;
forcing = zeros(400,1000);
for j = 1:400
    if max(df_Efin_F0(j,:))<0
        for jj=1:1000
            forcing(j,jj) = forcing_arctic+A*sin(jj*(2*3.141/1000)+B);
        end
    elseif min(df_Efin_F0(j,:))<0
        for jj=1:1000
            forcing(j,jj) = forcing_arctic+A*sin(jj*(2*3.141/1000)+B);
        end
    end
end

[tfin, ht_Efin_F0_arctic3,ht_Efinnt_F0_arctic3,ht_Tfin_F0_arctic3,ht_Tfinnt_F0_arctic3,ht_ffin_F0_arctic3,ht_ffinnt_F0_arctic3,ht_hfin_F0_arctic3,ht_hfinnt_F0_arctic3] = sea_ice_EBM_EEW22(3,0,Fb,dur,nt,n,forcing,df_Efinnt_F3);

%% Nudging method SEA ICE LOSS
% Only nudge areas where sea ice is present in control but not in future
nudge_tau = 7;
[tfin, ng_Efin_F0_nudge3,ng_Efinnt_F0_nudge3,ng_Tfin_F0_nudge3,ng_Tfinnt_F0_nudge3,ng_ffin_F0_nudge3,ng_ffinnt_F0_nudge3,ng_hfin_F0_nudge3,ng_hfinnt_F0_nudge3] = sea_ice_EBM_EEW22(4,0.0,Fb,dur+100,nt,n,df_Efinnt_F3,nudge_tau);
