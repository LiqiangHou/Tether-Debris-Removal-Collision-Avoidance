function plot_trajectory_and_thrust_profile

clear all;
close all;
global self
global MATRIX_A

format long g;


MATRIX_A = eye(7 , 7)*1.0e-12;



% constant paramters
dayszero = 0.15;    % p=0.10;
initialize_parameter(dayszero);

p0       = diag( [0.001 0.001 0.001 0.001 0.001 0.001 0.001]*1.0e-0);

self.p0  = p0;
self.tau = [0.0 , 1.0];


% the solution
% feq = [     0.00518461498391793      -0.00136506947682165       0.00746618131565968       -0.0100697825790903         0.244086808814658         0.960638646307336       0.441788421772173]
xoptm = [0.0383921584146486	0.00411304099871005	0.223911689663605	-0.448626642563718	-0.00233822057439055	0.441811560481595	0.0375187731455364];


% feq = [    0.00910619305097955      -0.00656041928460005
% 0       0.00544582565864205         0.692974539167654
% 0.912832126423306          0.29898300627914];
% ode113 option:RfEvents
xoptm = [0.0990177915105663       0.00917255395250577        0.0887312321784412        -0.103009615889113      -0.00614950306603373         0.298989623495269         0.118144838982202];
% settings of ode
tspan = [0, 1];

% 
tau_0 = tspan(1);
tau_f = tspan(2);
plot_fitness(tau_0,tau_f,xoptm);



return
end

%


function plot_fitness(tau_0,tau_f,x0)
global r0bar u0bar v0bar theta0 alfa0 w0bar ufbar alfaf wfbar lambdaf_alfa lambdaf_wbar l mu kappa lambdaf_theta
global rfbar
global tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0

global fitflux fitfgrad
global kappa_alpha alpha_threshold
global m1

global days
global max_delta_alpha

%State vecotr X has components
%x(1)=rbar, nondimensionalized radius from center of earth
%x(2)=theta
%x(3)=alfa
%x(4)=ubar
%x(5)=vbar
%x(6)=wbar
%x(7)=lambda_rbar
%x(8)=lambda_alfa
%x(9)=lambda_ubar
%x(10)=lambda_vbar
%x(11)=lambda_wbar
%x(12)=mass
%x(13)=lambda_mass

% initial value of state and costate of the ode, 
% x0: optimal initial value of cosctate
y0   = [r0bar theta0 alfa0 u0bar v0bar w0bar x0(1) x0(2) x0(3) x0(4) x0(5) m1 x0(6)];
days = x0(7);

% Initial amount of days to thrust. 
tf =24*3600*days; % Days in seconds


%
tspan   = [0 , 1] ;



atol = 1.0e-6;
rtol = 1.0e-6;


% options = odeset('RelTol',rtol,'AbsTol',atol,'Events',[]);
options = odeset('RelTol',rtol,'AbsTol',atol,'Events',@RfEvents);
%              
sol  = ode113(@odes_function_flux, tspan, y0,options);



Nt      = 10000;
tau_end = sol.xe;
tau     = linspace(0,tau_end,Nt);

plot_results(sol,tau);


return
end





function initialize_parameter(days)
global v0 u0 l l1 tf m12 m m2 p m1 r0 mu r0bar u0bar v0bar theta0 alfa0 w0bar ufbar lambdaf_ubar ...
lambdaf_alfa alfaf lambdaf_wbar wfbar
global rfbar
global fitflux  fitfgrad
global kappa
global area_c
global tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0
global ISP
global kappa_alpha alpha_threshold
global period
mu=3.98e5;
r0=7200;

u0=0;
v0    = (mu/r0)^0.5/r0;
alfa0 = pi/2 -deg2rad(5);


l=10;
t0=0;
theta0=0;
w0=0;

m1=1500;
m2=100;
m=m1+m2;
m12=m1*m2/m;
l1= l*m2/m;


% p=0.025;    

p=0.10;       % test case of 0.1N thrust 

ISP = 3000;

% Initial amount of days to thrust. 
tf =24*3600*days; % Days in seconds



% initial state
r0bar  = 1;
u0bar  = 0;
v0bar  = 1;
w0bar  = 0;
theta0 = 0;

% final conditions to be met
ufbar        =0;
lambdaf_alfa = 0;
lambdaf_wbar = 0;


rfbar = (6378+350)/r0;

period = 2*pi*r0/v0;

% debris flux data
SMA_data          = dlmread('Flux_altitude_data.txt');
x  = SMA_data(:,1); %SMA
y  = SMA_data(:,19); %flux


fitflux = fit(x,y,'cubicinterp');
Gradient_SMA_data=gradient(y,x);
fitfgrad = fit(x,Gradient_SMA_data,'cubicinterp');

% slack constraint constant 
kappa = 0.50*1.0e3;

% area of tether system,m^2
area_c=36*pi;

% probability constraint 
v_threshold = 1.0e-4;

% slack variable of alpha constraint
kappa_alpha     = 5;
alpha_threshold = deg2rad(5.0);




return
end







function J = Obj_Convolution(x)



% initial value of costate and days
xzero  = x(1:7);                


% initial variance of the final condition

tau_0 = 0.0;
tau_f = 1.0;      % tau_f should be keep unvaried
[feq]  = fitness(tau_0,tau_f,xzero);
Mu    = feq;



% SIGMA =  Jacobian_Covariance(tau_0,tau_f,xzero);
% 
          
% SIGMA: covraiance of initila guess  uisng Jacobian 
% 
% % % p=0.1N, sucess,deg45
SIGMA =[...
        0.0335873890035292        0.0123059851447506       -0.0224488706419929       0.00435895276899434        0.0255799397700252        0.0255799416027115        0.0528380895221596
        0.0123059851447506       0.00450875387684673      -0.00822497600536212       0.00159706394612208       0.00937215925835937       0.00937215992983193         0.019359192959845
       -0.0224488706419929      -0.00822497600536212        0.0150041967551564      -0.00291340201631715       -0.0170969157163996       -0.0170969169413158       -0.0353155000088981
       0.00435895276899434       0.00159706394612208      -0.00291340201631715      0.000565702479592174       0.00331975043608015       0.00331975067392516       0.00685729803548551
        0.0255799397700252       0.00937215925835937       -0.0170969157163996       0.00331975043608015        0.0194815178568766        0.0194815192526387        0.0402411496588206
        0.0255799416027115       0.00937215992983193       -0.0170969169413158       0.00331975067392516        0.0194815192526387        0.0194815206484009        0.0402411525419158
        0.0528380895221596         0.019359192959845       -0.0353155000088981       0.00685729803548551        0.0402411496588206        0.0402411525419158        0.0831223797735037
        ];
% % % p=0.1N, sucess!!!
s =         [0.01 0.01 0.001 0.001 0.01 0.0005 0.001]*5.0e1   ;    ...  %        variance of boundary  sucess !!!



S = diag(s);

% performance index
M = zeros(1,7);
J = quadratic_convolution(Mu,M,SIGMA,S);


return
end








function [value,isterminal,direction] = RfEvents(t,y)

global rfbar
value = y(1) - rfbar;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;     % The zero can be approached from either direction
return
end




function plot_results(sol,tau)

global v0 u0 l l1 tf m12 m m2 p m1 r0 mu r0bar u0bar v0bar theta0 alfa0 w0bar ufbar lambdaf_ubar ...
lambdaf_alfa alfaf lambdaf_wbar wfbar

global period


% result structure of the ODE
Z     = deval(sol,tau);


r_sol     = Z(1,:)*r0; % Radius [km]
theta_sol = Z(2,:); % Radial component of velocity [km/s]
alfa_sol  = Z(3,:); % Tangential component of velocity [km/s]
u_sol     = Z(4,:)*v0*r0;
v_sol     = Z(5,:)*v0;
w_sol     = Z(6,:)*v0;
lambda_rbar_sol = Z(7,:);
lambda_alfa_sol = Z(8,:);
lambda_ubar_sol = Z(9,:);
lambda_vbar_sol = Z(10,:);
lambda_wbar_sol = Z(11,:);

%%plot and results
final_radius = r_sol(end)
final_theta = theta_sol(end);
X = r_sol.*cos(theta_sol);
Y = r_sol.*sin(theta_sol);
time = tau*(tf);

beta = atand((lambda_ubar_sol./r0 +...
lambda_wbar_sol.*m2.*sin(alfa_sol)./(m12.*l))./(lambda_wbar_sol.*m2.*cos(alfa_sol)./(m12*l)...
    +lambda_wbar_sol./(r0.*r_sol)-lambda_vbar_sol./(r0.*r_sol)));
A = alfa_sol;
beta1 = asind ((lambda_ubar_sol./r0 ...
    +lambda_wbar_sol.*m2.*sin(alfa_sol)./(m12.*l))./(((lambda_ubar_sol./r0 ...
    +lambda_wbar_sol.*m2.*sin(alfa_sol)./(m12.*l)).^2+(lambda_wbar_sol.*m2.*cos(alfa_sol)./(m12*l)...
+lambda_wbar_sol./(r0.*r_sol)-lambda_vbar_sol./(r0.*r_sol)).^2).^0.5));

figure(1)
plot(X,Y,'k')
xlabel('x-direction [km]','fontsize',14)
ylabel('y-direction [km]','fontsize',14)
title('Trajectory of the Tethered System','fontsize',16)
axis equal
hold on
plot(0,0,'k.','MarkerSize',20)
plot(r0*cos(theta0),r0*sin(theta0),'-.rx','MarkerSize',10)
plot(final_radius*cos(final_theta),final_radius*sin(final_theta),'-.ro','MarkerSize',10)
% figure(3)
% plot(X1,Y1)
figure(2)
plot(time/period,beta1,'k')
xlabel('time [number of inital period]','fontsize',14)
ylabel('Thrust Angle [deg]','fontsize',14)
title('Thrust Angle Vs Time','fontsize',16)
figure(3)
plot(time/period,A*180/pi,'k')
xlabel('time [number of inital period]','fontsize',14)
ylabel('Tether Angle [deg]','fontsize',14)
title('Tether Angle Vs Time','fontsize',16)
figure(4)
plot(time/period,beta,'k')
xlabel('time [number of inital period]','fontsize',14)
ylabel('Thrust Angle [deg]','fontsize',14)
title('Thrust Angle Vs Time','fontsize',16)
tf

figure(5)
alt = r_sol - 6378.0;

plot(time/period,alt,'k')
xlabel('time [number of inital period]','fontsize',14)
ylabel('Altitude [km]','fontsize',14)
title('Altitude Vs Time','fontsize',16)
tf





return
end



