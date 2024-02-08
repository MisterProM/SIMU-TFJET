function y=rev_h(T,h0,alpha)
%With this function, we can obtain the properties of the mass flux for a
%specific alpha and temperature.

r=287.15; %J/KgK
%Entalpy
y=r*((3.5*T-1.4*10^-5*T^2+7.467*10^-9*T^3+3090/(exp(3090/T)-1))+alpha*(-149.054+4.47659*T+4.00997*10^-3*T^2-6.12432*10^-7*T^3))/(1+alpha)-h0;

end
