function [flux,f_grad]=flux_valueandpartial_2d(x0,fitflux, fitfgrad) 




x0 = x0*7.171*10^3 ;%convert r bar to radius in kilometer

%locate the point
flux = fitflux(x0);
flux = flux/365;% convert flux from 1/m^2/yr to 1/m^2/day

%gradient portion
f_grad = fitfgrad(x0);
% f_grad = f_grad/365*1000; % convert flux partial from 1/m^2/yr to 1/m^2/day also makes 
% %partial per m of altitude rather than per km of altitude

f_grad = f_grad/365; % convert flux partial from 1/m^2/yr to 1/m^2/day also makes 


%return the flux  and flux partial at a radius x0
end
