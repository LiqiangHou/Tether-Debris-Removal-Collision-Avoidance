function dX_dtau=odes_function_flux(tau,x)
%Independent variable here is nondimensional time, tao, and the state
%vector is X
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

global tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0
global fitflux  fitfgrad
global kappa
global area_c
global ISP
global max_delta_alpha
global kappa_alpha


m1 = x(12);
m=m1+m2;
m12=m1*m2/m;
l1= l*m2/m;


% % % flux and its gradient
% [flux,flux_partial]=flux_valueandpartial_2d(x(1),fitflux, fitfgrad) ;
% 



% %

% beta = atan((x(9)/r0 + x(11)*m2*sin(x(3))/(m12*l))/(x(11)*m2*cos(x(3))/(m12*l)+x(11)/(r0*x(1))-x(10)/(r0*x(1))));

denominator = x(11)*m2*cos(x(3))/(m12*l)+x(11)/(r0*x(1))-x(10)/(r0*x(1));
if(abs(denominator) < 1.0e-2)
    beta_denominator = sign(denominator)*1.0e-2;
else
    beta_denominator = denominator;
end
beta = atan((x(9)/r0 + x(11)*m2*sin(x(3))/(m12*l))/(beta_denominator));

rho = switch_function(x,beta);
if(rho > 0)
    const_p = p;
    p = 0;
else
    p = p;
end

drbar_dtau=tf*v0*x(4);
dtheta_dtau=tf*v0*x(5);
dalfa_dtau=tf*v0*x(6);


dubar_dtau=tf*p*sin(beta)/(v0*r0*m)+tf*x(1)*v0*x(5)^2-mu*tf/(v0*r0^3*x(1)^2)... 
    -3*mu*m12*l^2*tf*(1-3*(cos(x(3)))^2)/(2*m*r0^5*v0*x(1)^4);
dvbar_dtau = 3*mu*m12*tf*l^2*sin(2*x(3))/(2*m*v0*r0^5*x(1)^5) - 2*v0*x(4)*x(5)*tf/x(1)...
    -cos(beta)*p*tf/(m*r0*v0*x(1));
dwbar_dtau = m2*cos(beta-x(3))*p*tf/(v0*m*m12*l) - 3*mu*tf*sin(2*x(3))/(2*v0*r0^3*x(1)^3)...
    -3*mu*m12*tf*l^2*sin(2*x(3))/(2*m*v0*r0^5*x(1)^5) + 2*v0*x(4)*x(5)*tf/x(1)... 
+cos(beta)*p*tf/(m*r0*v0*x(1));

dmass_datu  = -tf*p/ISP;


%% dlambda_dtau

lambda_theta = 0.0;

dlambda_rbar_dtau =...
x(11)*((2*tf*x(4)*v0*x(5))/x(1)^2 - (9*mu*tf*sin(2*x(3)))/(2*r0^3*x(1)^4*v0) + (p*tf*cos(beta))/(m*r0*x(1)^2*v0) - (15*l^2*m12*mu*tf*sin(2*x(3)))/(2*m*r0^5*x(1)^6*v0)) - x(10)*((2*tf*x(4)*v0*x(5))/x(1)^2 + (p*tf*cos(beta))/(m*r0*x(1)^2*v0) - (15*l^2*m12*mu*tf*sin(2*x(3)))/(2*m*r0^5*x(1)^6*v0)) - x(9)*(tf*v0*x(5)^2 + (2*mu*tf)/(r0^3*x(1)^3*v0) - (6*l^2*m12*mu*tf*(3*cos(x(3))^2 - 1))/(m*r0^5*x(1)^5*v0));
 
 
dlambda_ubar_dtau =... 
(2*x(10)*tf*v0*x(5))/x(1) - x(7)*tf*v0 - (2*x(11)*tf*v0*x(5))/x(1);
 
 
dlambda_vbar_dtau =... 
(2*x(10)*tf*x(4)*v0)/x(1) - 2*x(9)*x(1)*tf*v0*x(5) - lambda_theta*tf*v0 - (2*x(11)*tf*x(4)*v0)/x(1);
 
 

% % p=0.1N, sucess,deg45

% delta_dlambda_alfa_dtau due to constraing of alpha
if(abs(x(3) - pi/2) <= deg2rad(45))
    delta_dlambda_alfa_dtau = 0.0;
else
    if((x(3) - pi/2) >= deg2rad(45))
        delta_dlambda_alfa_dtau = kappa_alpha;
    end
    if((x(3) - pi/2) < -deg2rad(45))
        delta_dlambda_alfa_dtau = -kappa_alpha;
    end
end



% dlambda_alfa_dtau =... 
% x(11)*((3*mu*tf*cos(2*x(3)))/(r0^3*x(1)^3*v0) + (m2*p*tf*sin(x(3) - beta))/(l*m*m12*v0) + (3*l^2*m12*mu*tf*cos(2*x(3)))/(m*r0^5*x(1)^5*v0)) - (3*l^2*x(10)*m12*mu*tf*cos(2*x(3)))/(m*r0^5*x(1)^5*v0) + (9*l^2*x(9)*m12*mu*tf*cos(x(3))*sin(x(3)))/(m*r0^5*x(1)^4*v0);

dlambda_alfa_dtau =... 
x(11)*((3*mu*tf*cos(2*x(3)))/(r0^3*x(1)^3*v0) + (m2*p*tf*sin(x(3) - beta))/(l*m*m12*v0) + (3*l^2*m12*mu*tf*cos(2*x(3)))/(m*r0^5*x(1)^5*v0)) - (3*l^2*x(10)*m12*mu*tf*cos(2*x(3)))/(m*r0^5*x(1)^5*v0) + (9*l^2*x(9)*m12*mu*tf*cos(x(3))*sin(x(3)))/(m*r0^5*x(1)^4*v0)...
+ delta_dlambda_alfa_dtau;
 
 
dlambda_wbar_dtau = -x(8)*tf*v0;
 
 
dlambda_mass_dtau =... 
x(11)*((p*tf*cos(beta))/(m^2*r0*x(1)*v0) + (m2*p*tf*cos(x(3) - beta))/(l*m^2*m12*v0) - (3*l^2*m12*mu*tf*sin(2*x(3)))/(2*m^2*r0^5*x(1)^5*v0)) + x(9)*((p*tf*sin(beta))/(m^2*r0*v0) + (3*l^2*m12*mu*tf*(3*cos(x(3))^2 - 1))/(2*m^2*r0^5*x(1)^4*v0)) - x(10)*((p*tf*cos(beta))/(m^2*r0*x(1)*v0) - (3*l^2*m12*mu*tf*sin(2*x(3)))/(2*m^2*r0^5*x(1)^5*v0));



%%
dX_dtau = [drbar_dtau; dtheta_dtau; dalfa_dtau; dubar_dtau; dvbar_dtau; dwbar_dtau;dlambda_rbar_dtau;dlambda_alfa_dtau;dlambda_ubar_dtau; dlambda_vbar_dtau; dlambda_wbar_dtau; dmass_datu; dlambda_mass_dtau;];


% 
% 
if(rho > 0)
    % restore the constant value of p
    p = const_p;
else
    p = p;
end
return
end



function rho = switch_function(x,beta)
global flux flux_partial tf_days t0_days area_c v_threshold delta_t tf v0 p mu m12 l m m2 r0
global ISP


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



% rho= (lambda_ubar / r_0) (sin(beta) / m)  ...
%     - lambda_vbar * cos(beta)/ (m * r0 * rbar)} ...
%     + lambda_wbar * ((m2 * cos(beta-alpha))/(m * m12 * l) + cos(beta) / (m * r0 * rbar))
%     - lambda_mass * (v0 /Isp)



rho= x(9)/r0*sin(beta)/m             ...
   - x(10)* cos(beta)/(m*r0*x(1))  ...
   + x(11)* (m2*cos(beta-x(3))/(m*m12*l) + cos(beta)/(m*r0*x(1))) ...
   - x(13) * (v0 /ISP);








return
end
 