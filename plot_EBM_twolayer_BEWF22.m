% This code plots output from the code two_ocean_layer_EBM_2022.m. If using
% the default parameter values given below, the code produces a plot 
% similar to the top row of Figure 3 in Beer, Eisenman, Wagner and Fine 
% (2022; see reference)
%
% E. Beer, I. Eisenman, T.J.W. Wagner and E.C. Fine (2022). A possible 
% hysteresis in the Arctic Ocean due to release of subsurface heat during 
% sea ice retreat. Submitted to Journal of Physical Oceanography.

% Default parameter case
clearvars
F_in = zeros(1,2000); % run F = 0 radiative forcing for 2000 years (W m^-2)
kv_w = 5;             % vertical heat flux coefficient when T > 0 (W m^-2 K^-1)
kv_i = 2;             % vertical heat flux coefficient when T < 0 (W m^-2 K^-1)
Ds = 0.5;             % diffusion coefficient for heat transport in mixed layer (W m^-2 K^-1)
Dd = 0.15;            % diffusion coefficient for heat transport in deep layer (W m^-2 K^-1)
A = 192;              % OLR when T = 0 (W m^-2)


% Run model
[x_out,~,T_out,Td_out,Fb_out] = two_ocean_layer_EBM_2022(F_in,kv_w,kv_i,Ds,Dd,A);
disp('default run completed')

% Take the last value
T_end = T_out(end,:);
Td_end = Td_out(end,:);
Fb_end = Fb_out(end,:);

% Calculate the surface mixed layer temperature
T_sml = T_end;
T_sml(T_sml<-2) = -2;

% Set up the figure 
f1 = figure; 
f1.Position = [54 361 1172 405];

% Plot the surface mixed layer temperature
subplot(1,3,1), hold on
plot(asind(x_out),T_sml,'.-')
xlim([0 90])
xlabel('Latitude');
ylabel('SML temperature (\circC)')
set(gca,'fontsize',14)

% Plot the deep layer temperature
subplot(1,3,2), hold on
plot(asind(x_out),Td_end,'.-')
xlim([0 90])
xlabel('Latitude');
ylabel('DL temperature (\circC)')
set(gca,'fontsize',14)

% Plot the  vertical heat flux temperature
subplot(1,3,3),hold on
plot(asind(x_out),Fb_end,'.-')
xlim([0 90])
xlabel('Latitude');
ylabel('Vertical heat flux (W m^{-2})')
set(gca,'fontsize',14)