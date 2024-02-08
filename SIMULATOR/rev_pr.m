function y= rev_pr(T,phi,alpha)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Tr=273.15;
%phi_airTr=3.5*log(Tr)-2.8*10^-5*Tr+1.12*10^-8*Tr^2+3090/(Tr*(exp(3090/Tr)-1))-log((exp(3090/Tr)-1)/exp(3090/Tr));
%phi_air=3.5*log(T)-2.8*10^-5*T+1.12*10^-8*T^2+3090/(T*(exp(3090/T)-1))-log((exp(3090/T)-1)/exp(3090/T));
%phi_fuel=4.47659*log(T)+8.01994*10^-3*T+9.18648*10^-7*T^2;
%y=(phi_air+alpha*phi_fuel)/(1+alpha)-(phi_airTr+alpha*phi_fuel)/(1+alpha)-phi;
y=(3.5*log(T)-2.8*10^-5*T+1.12*10^-8*T^2+3090/(T*(exp(3090/T)-1))-log((exp(3090/T)-1)/exp(3090/T))+alpha*(4.47659*log(T)+8.01994*10^-3*T+9.18648*10^-7*T^2))/(1+alpha)-phi;

end