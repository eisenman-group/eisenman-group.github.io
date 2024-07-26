% This code computes the differential and effective net feedback
% parameters given the surface temperature (T) and the
% net radiative response of the climate system (F_{net} = N - F{GHG}),
% with all quantities averaged annually and globally.
%
% This code is shared with the following paper:
% I. Eisenman and K. Armour (2024). The radiative feedback continuum from
% Snowball Earth to an ice-free hothouse. Nature Communications,
% https://doi.org/10.1038/s41467-024-50406-w.
%
% The code generates the values plotted in Figs. 2d,e of the paper given
% the simulation output for F_{net}(t) and T(t) (plotted in Fig. 2c). For
% F_{net} and T, it requires the source data for the line plots in the
% paper, which is provided as a supplemental data zip-file with the paper
% and can also be downloaded from
%   https://eisenman.ucsd.edu/code/source_data_EA24.zip
% and extracted. This code reads in the file
%   data/fig2.txt
%
% Ian Eisenman (eisenman@ucsd.edu), 2024

% Read in time series data for T(t) and Fnet(t)
d = load('data/fig2.txt');
T=d(:,1);
Fnet=d(:,3)-d(:,2);
tw=(1:279)'; % time during Warming run
tc=(1:514)'; % time during Cooling run
t=[-flipud( tc ); tw]; % single time array corresponding to T and Fnet

% surface temperature and radiative response averaged over years 480-499 of the PI simulation
T0 = 15.0523; F0 = 0.0543;

% effective net feedback
p=constrained_polyfit_fn(T-T0,Fnet-F0,12,0,0);
lambda_eff=polyval(p,T-T0);

% differential net feedback
lambda_diff=nan(size(T,1),1);
for m=1:length(T)
    g1=find(T<T(m)-3,1,'last');
    g2=find(T>T(m)+3,1,'first');
    if ~isempty(g1) && ~isempty(g2)
        p2 = tls_regression_fn(Fnet(g1:g2),T(g1:g2),t(g1:g2));
        lambda_diff(m) = p2(1);
    end
end

figure(1), clf
[~,g]=sort(T);
subplot(1,2,1)
plot(T(g),lambda_diff(g))
xlabel('T (^oC)'), ylabel('\lambda_{net}^{diff} (W/m^2/K)')
subplot(1,2,2)
plot(T(g),lambda_eff(g))
xlabel('T (^oC)'), ylabel('\lambda_{net}^{eff} (W/m^2/K)')

function p = constrained_polyfit_fn(x,y,n,x0,y0)
% polyfit(x,y,n) but constrained to go through points (x0,y0).
% Code adapted from
% https://www.mathworks.com/matlabcentral/answers/94272-how-do-i-constrain-a-fitted-curve-through-specific-points-like-the-origin-in-matlab

x = x(:); %reshape the data into a column vector
y = y(:);
x0=x0(:); y0=y0(:);
% 'C' is the Vandermonde matrix for 'x'
% 'n' is the degree of polynomial to fit
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end
C = V;
% 'd' is the vector of target values, 'y'.
d = y;
% There are no inequality constraints in this case, i.e.,
A = [];
b = [];
% We use linear equality constraints to force the curve to hit the required point. In
% this case, 'Aeq' is the Vandermoonde matrix for 'x0'
Aeq = x0.^(n:-1:0);
% and 'beq' is the value the curve should take at that point
beq = y0;
p = lsqlin( C, d, A, b, Aeq, beq );
p = p(1:end-1);
end

function b = tls_regression_fn(y,x,t)
% Total least squares estimate of the regression of y on x [linear
% relationship between x(t) and y(t)]. This is the maximum
% likelihood estimator of the relationship when the units of x and y are
% chosen to be the standard deviations of the respective
% normally-distributed error terms.
% Here the error terms are estimated using the residuals of the time
% trends. This is equivalent to eqs. (2) & (3) in
% Winton (2010, https://doi.org/10.1175/2011JCLI4146.1).
% Code by Ian Eisenman, 2010, updated 2013 and 2024

dx=x-mean(x);
dy=y-mean(y);
n=length(dx)-1;

sx=sqrt(sum(dx.*dx)/n); % std(x)
sy=sqrt(sum(dy.*dy)/n); % std(y)
rhoyx=sum(dy.*dx)/(n*sx*sy);

dt=t-mean(t);
st=sqrt(sum(dt.*dt)/n); % std(t)
rhoxt=sum(dx.*dt)/(n*sx*st);
rhoyt=sum(dy.*dt)/(n*sy*st);
L = (sy^2/sx^2) * (1-rhoyt^2)/(1-rhoxt^2);

b = ( sy^2-L*sx^2+sqrt( (sy^2-L*sx^2)^2 + 4*L*rhoyx^2*sy^2*sx^2 ) ) / ( 2*rhoyx*sy*sx ) ;
end
