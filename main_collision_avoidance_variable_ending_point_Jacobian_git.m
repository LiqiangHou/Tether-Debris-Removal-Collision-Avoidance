function main_tether_debris_removal_collision_avoidance_git

clear all;
close all;
global self
global MATRIX_A
global SIGMA

format long g;



MATRIX_A = eye(7 , 7)*1.0e-13;


% constant paramters
dayszero = 0.12;    
initialize_parameter(dayszero);

p0       = diag( [0.001 0.001 0.001 0.001 0.001 0.001 0.001]*1.0e-0);

self.p0  = p0;
self.tau = [0.0 , 1.0];

% bounds of expected vlaue of the initila costate
daysmin  = 0.02;
daysmax  = 0.30;    

lb      = [ones(1,6) * -10   daysmin];
ub      = [ones(1,6) *  10   daysmax];


% the second part of design variables: varaiance of the terminal boundary
s_lb = 1.0e-8 * [ones(1,7)];
s_ub = 5.0e2  * [ones(1,7)];

s_0 = [5.0e-1*ones(1,4) , 0.10e-1, 0.10e-1, 0.10e-1];     



% objective function with convolution
fun = @Obj_Convolution;

A   = [];
b   = [];
Aeq = [];
beq = [];
nonlcon = [];

% Jocabian of the guesses  
x0     = [ones(1,6)*0.1 , dayszero];
SIGMA  =  Jacobian_Covariance(0,1,x0);



% Design variables: Initial guesses and varaiance of the terminal boundary

x_s_0 =  [x0 , s_0];
x_s_lb = [lb , s_lb];
x_s_ub = [ub , s_ub];

   
options      = optimoptions('fmincon','Algorithm','sqp','Display','iter',...
                        'StepTolerance',1.0e-30,'OptimalityTolerance',1.0e-30,...
                        'MaxFunctionEvaluations',20000,'OutputFcn',@outfun);   

[xoptm,fval] = fmincon(fun,x_s_0,A,b,Aeq,beq,x_s_lb,x_s_ub,nonlcon,options);



% x0: optimal initial value of cosctate
% x0 = [     -0.000843663892999919     -0.000843663892999919     -0.000843663892999919     -0.000843663892999919     -0.000843663892999919            0.196889010407           0.1474294752341];


return
end

%


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


% test case of 0.1N thrust 
p=0.10;       

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








%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function stop = outfun(x,optimValues,state) 
global self
global SIGMA
stop   = false;


tau    = self.tau;
tau_0  = tau(1);
tau_f  = tau(end);

xzero  = x(1:7)';
[feq]  = fitness(tau_0,tau_f,xzero)

if norm(feq) < 0.01 
    stop = true;
    disp('Stopping, norm(ceq) < 0.1')
end

SIGMA =  Jacobian_Covariance(tau_0,tau_f,xzero);


return
end


function J = Obj_Convolution(x)

global SIGMA
global s

% initial value of costate and days
xzero  = x(1:7);                
s      = x(8:14);

% initial variance of the final condition

tau_0 = 0.0;
tau_f = 1.0;      % tau_f should be keep unvaried
[feq]  = fitness(tau_0,tau_f,xzero);
Mu    = feq;






S = diag(s);

% performance index
M = zeros(1,7);
J = quadratic_convolution(Mu,M,SIGMA,S);


return
end







function ceq = fitness(tau_0,tau_f,x0)
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
y0   = [r0bar theta0 alfa0 u0bar v0bar w0bar x0(1) x0(2) x0(3) x0(4) x0(5) m1 x0(6)];
days = x0(7);

% Initial amount of days to thrust. 
tf =24*3600*days; % Days in seconds


%
tspan   = [0 , 1] ;

atol = 1.0e-6;
rtol = 1.0e-6;


options = odeset('RelTol',rtol,'AbsTol',atol,'Events',@RfEvents);

%              
[t,Y]   = ode113(@odes_function_flux, tspan, y0,options);

Yf = Y(end,:);

% flux and probability
[flux,flux_partial]=flux_valueandpartial_2d(Yf(1),fitflux, fitfgrad) ;
prob    = 1 - exp(-flux * area_c * days);
probdot = 1 - exp(-flux * area_c);

dYf = odes_function_flux(tau_f,Yf); 

lambdaf  = [Yf(7:11),Yf(13)];
dxf      = [dYf(1) , dYf(3) , dYf(4) , dYf(5) , dYf(6) , dYf(12)]';      % theta is not considered

days_f   = t(end)*days; 
alphaf   = Yf(3);

% Hamilton
% constriant of alpha: [90-5deg , 90+5deg], constriant of collsion probability:  1.0e-4

Hf       = 1 + lambdaf*dxf + kappa*max(0 , (probdot - v_threshold/(days_f * 86400))) + kappa_alpha*max(0 , (abs(alphaf - pi/2) - deg2rad(5))); 



% Terminal boundary condition
ceq  = [
Yf(8)  - lambdaf_alfa
Yf(11) - lambdaf_wbar
Yf(1)  - rfbar             % Final Condition
Yf(5)  - sqrt(1/rfbar^3)   % Final Condition
-Yf(7) + 3/2*Yf(10)/rfbar^(5/2) + 1                                       %final
Hf                         % Final condition of minmum ToF, rfbar and vfbar are time unvaried
Yf(13)                     % lambda_massf
]';



if(Yf(1) < 0.90)
    disp(ceq);
end



return
end

function SIGMA =  Jacobian_Covariance(tau_0,tau_f,x0)
ceq_zero = fitness(tau_0,tau_f,x0);
x        = x0; 



deltax      = ones(length(x0),1)*1.0e-5;

% 



for i=1:length(x0);
    x(i)= x0(i) + deltax(i);
    ceq = fitness(tau_0,tau_f,x);

    F(i,:) = (ceq - ceq_zero)'/deltax(i);
end

SIGMA = (F*deltax)*(F*deltax)';

return
end


function [value,isterminal,direction] = RfEvents(t,y)

global rfbar
value = y(1) - rfbar;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;     % The zero can be approached from either direction
return
end

function  J = quadratic_convolution(Mu,M,SIGMA,S)
global MATRIX_A

% A = MATRIX_A;
% J = (Mu - M) * pinv(SIGMA + A, 1.0e-20) * (Mu - M)' + trace(pinv(SIGMA + A, 1.0e-20) * S);
% J = (Mu - M) * invChol_mex(SIGMA + A) * (Mu - M)' + trace(invChol_mex(SIGMA + A) * S);



% A = 1.0e-12*eye(size(SIGMA));
 A = MATRIX_A;
J = (Mu - M) * pinv(SIGMA + A, 1.0e-20) * (Mu - M)' + trace(pinv(SIGMA + A, 1.0e-20) * S);
% J = (Mu - M) * invChol_mex(SIGMA + A) * (Mu - M)' + trace(invChol_mex(SIGMA + A) * S);


return
end

