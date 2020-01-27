% This code plots output from the code sea_ice_EBM_deep_BEW20.m. If no 
% climate forcing is input (F=0), the code produces a plot similar to the 
% top row of Figure 2 in Beer, Eisenman, and Wagner (2020; see reference).
%
% E. Beer, I. Eisenman, and T.J.W. Wagner (2020). Polar amplification due 
% to enhanced heat flux across the halocline. Geophys Res Lett.
%
%--------------------------------------------------------------------------

%%Run model code
F = 0; % default climate forcing
[t, E0, ~, T0, Td0] = sea_ice_EBM_deep_BEW20(F); %output variables needed 

%%Define parameters
x = 1/200:1/100:1-1/200;
lat = asind(x); % latitude N
Lf = 9.5;     % sea ice latent heat of fusion (W yr m^-3)

% compute seasonal ice edge
xi = zeros(1,100);
for j = 1:length(t)
    Ej = E0(:,j);
    xe=x(find(Ej<0,1,'first'));
    if ~isempty(xe)
        xi(j) = xe;
    else
        xi(j) = max(x);
    end
end
lati = asind(xi);

%%Set up figure 
fig1 = figure; clf
fig1.Position = [237 371 955 394];

% color axis limits for temperature
colorvec = [-2 30];

% calculate surface mixed layer temperature (T_sml)
Tsml = T0;
Tsml(Tsml<-2) = -2;

% plot surface mixed layer temperature (Fig. 2a)
clevs = -2:4:40;
subplot(1,3,1)
contourf(t,lat,Tsml,clevs)
caxis(colorvec)
xlabel('Time (years)')
ylabel('Latitude')
ylim([0 lat(end)])
set(gca,'fontsize',14)
h = colorbar; set(get(h,'label'),'string','T_{sml} (\circC)');
h.Location = 'southoutside';
h.FontSize = 14;

% set colorbar for temperature subplots
map = parula(8);
colormap(gca,map)

% plot the deep layer temperature (Fig. 2b)
subplot(1,3,2)
contourf(t,lat,Td0,clevs) % use the same contours (clevs) as for T_sml
caxis(colorvec)
xlabel('Time (years)')
ylim([0 lat(end)])
set(gca,'fontsize',14)
h = colorbar; set(get(h,'label'),'string','T_{d} (\circC)');
h.Location = 'southoutside';
h.FontSize = 14;
colormap(gca,map)

% calculate the ice thickness using the surface enthalpy 
hfin = -(E0/Lf.*(E0<0));

% plot the ice thickness (Fig. 2c)
clevsice = 0.0001:.6:4; % set different contours for the ice subplot
subplot(1,3,3)
contourf(t,lat,hfin,clevsice)
caxis([0 3.6])
xlabel('Time (years)')
ylim([0 lat(end)])
set(gca,'fontsize',14)
h = colorbar; set(get(h,'label'),'string','h (m)');
h.Location = 'southoutside';
h.FontSize = 14;

% set colorbar for ice subplot
mapice = parula(6);
colormap(gca,mapice)

