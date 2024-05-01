function plot_trajectory_and_thrust_profile_git

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
xoptm = [0.0990177915105663       0.00917255395250577        0.0887312321784412        -0.103009615889113      -0.00614950306603373         0.298989623495269         0.118144838982202];


% settings of ode
tspan = [0, 1];
tau_0 = tspan(1);
tau_f = tspan(2);

%
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

%State vecotr components
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
%x(14)=collison_prob

% initial value of state and costate of the ode, 
% x0: optimal initial value of cosctate
y0   = [r0bar theta0 alfa0 u0bar v0bar w0bar x0(1) x0(2) x0(3) x0(4) x0(5) m1 x0(6) 0.0];
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
sol  = ode113(@odes_function_flux_prob, tspan, y0,options);



Nt      = 10000;
tau_end = sol.xe;
tau     = linspace(0,tau_end,Nt);

plot_results(sol,tau,days);


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










function [value,isterminal,direction] = RfEvents(t,y)

global rfbar
value = y(1) - rfbar;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;     % The zero can be approached from either direction
return
end




function plot_results(sol,tau,days)

global v0 u0 l l1 tf m12 m m2 p m1 r0 mu r0bar u0bar v0bar theta0 alfa0 w0bar ufbar lambdaf_ubar ...
lambdaf_alfa alfaf lambdaf_wbar wfbar

global period
global fitflux  fitfgrad

global area_c



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
collision_prob  = Z(14,:);

%%plot and results
final_radius = r_sol(end)
final_theta = theta_sol(end);
X = r_sol.*cos(theta_sol);
Y = r_sol.*sin(theta_sol);
time = tau*(tf);       % tf equals days*86400

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
plot(time,beta1,'k')
xlabel('time [sec]','fontsize',14)
ylabel('Thrust Angle [deg]','fontsize',14)
title('Thrust Angle Vs Time','fontsize',16)
figure(3)
plot(time,A*180/pi,'k')
xlabel('time [sec]','fontsize',14)
ylabel('Tether Angle [deg]','fontsize',14)
title('Tether Angle Vs Time','fontsize',16)
figure(4)
plot(time,beta,'k')
xlabel('time [sec]','fontsize',14)
ylabel('Thrust Angle [deg]','fontsize',14)
title('Thrust Angle Vs Time','fontsize',16)
tf

figure(4)
alt = r_sol - 6378.0;

plot(time,alt,'k')
xlabel('time [sec]','fontsize',14)
ylabel('Altitude [km]','fontsize',14)
title('Altitude Vs Time','fontsize',16)


% flux data


r_data = Z(1,:);

for i = length(r_data)

    [flux,flux_partial]=flux_valueandpartial_2d(r_data(i),fitflux, fitfgrad) ;
    flux_data(i) = flux;

end

% plot collison
figure(5)
plot(time,collision_prob,'k')
xlabel('time [sec]','fontsize',14)
ylabel('','fontsize',14)
title('Accumulated Collision Probablity Vs Time','fontsize',16)



beta_set_rad = deg2rad(beta);
x_set        = Z;
[p_history,rho_histoy] = Thrust_profile(time,x_set,beta_set_rad);


figure(6)
% Thrust vs. time


plot(time,p_history,'k')
xlabel('time [sec]','fontsize',14)
ylabel('Thrust(N) ','fontsize',14)
title('Thrust Vs Time','fontsize',16)
ylim([-0.001,0.11]);





return
end


function  [p_history,rho_histoy] = Thrust_profile(time,x_set,beta_set_rad)
global tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0



for i=1:length(time)
    x        = x_set(:,i);
    beta_rad = beta_set_rad(i);

    rho = switch_function(x,beta_rad);
    if(rho > 0)
        const_p = p;
        p = 0;
    else
        p = p;
    end

    p_history(i)  = p;
    rho_histoy(i) = rho;


    % 
    % 
    if(rho > 0)
        % restore the constant value of p
        p = const_p;
    else
        p = p;
    end
end
return
end

function rho = switch_function(x,beta)
global flux flux_partial tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0
global ISP
global v0 u0 l l1 tf m12 m m2 p m1 r0 mu r0bar u0bar v0bar theta0 alfa0 w0bar ufbar lambdaf_ubar ...
lambdaf_alfa alfaf lambdaf_wbar wfbar


rbar = x(1); %nondimensionalized radius from center of earth
theta= x(2) ;
alfa= x(3)  ;
ubar= x(4) ;
vbar= x(5) ;
wbar= x(6) ;
lambda_rbar = x(7) ;
lambda_alfa = x(8);
lambda_ubar = x(9);
lambda_vbar= x(10);
lambda_wbar= x(11);
mass       = x(12) ;
lambda_mass= x(13) ;



% rho= (lambda_ubar / r_0) *(sin(beta) / m)  ...
%     - lambda_vbar * cos(beta)/ (m * r0 * rbar) ...
%     + lambda_wbar * ((m2 * cos(beta-alpha))/(m * m12 * l) + cos(beta) / (m * r0 * rbar))
%     - lambda_mass * (v0 /Isp);



rho= x(9)/r0*sin(beta)/m             ...
   - x(10)* cos(beta)/(m*r0*x(1))  ...
   + x(11)* (m2*cos(beta-x(3))/(m*m12*l) + cos(beta)/(m*r0*x(1))) ...
   - x(13) * (v0 /ISP);








return
end
 

