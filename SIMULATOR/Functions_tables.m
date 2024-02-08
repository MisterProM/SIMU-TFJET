function[h,phi,Cp,hf]=Functions_tables(T,alpha)
%With this function, we can obtain the properties of the mass flux for a
%specific alpha and temperature.

%Constant pressure specific heat
r=287.15; %J/KgK
Cp_air=3.5-2.8*10^-5*T+2.24*10^-8*T^2+(3090/T)^2*(exp(3090/T)/(exp(3090/T)-1)^2);
Cp_fuel=4.47659+8.01994*10^-3*T-1.8373*10^-6*T^2;
Cp=r*(Cp_air+alpha*Cp_fuel)/(1+alpha);

%Entalpy
h_air=3.5*T-1.4*10^-5*T^2+7.467*10^-9*T^3+3090/(exp(3090/T)-1);
h_fuel=-149.054+4.47659*T+4.00997*10^-3*T^2-6.12432*10^-7*T^3;
h=r*(h_air+alpha*h_fuel)/(1+alpha);

%Pressure ratio
Tr=273.15;
phi_airTr=3.5*log(Tr)-2.8*10^-5*Tr+1.12*10^-8*Tr^2+3090/(Tr*(exp(3090/Tr)-1))-log((exp(3090/Tr)-1)/exp(3090/Tr));
phi_air=3.5*log(T)-2.8*10^-5*T+1.12*10^-8*T^2+3090/(T*(exp(3090/T)-1))-log((exp(3090/T)-1)/exp(3090/T));
phi_fuel=4.47659*log(T)+8.01994*10^-3*T+9.18648*10^-7*T^2;
%phi=(phi_air+alpha*phi_fuel)/(1+alpha)-(phi_airTr+alpha*phi_fuel)/(1+alpha);
phi=(phi_air+alpha*phi_fuel)/(1+alpha);
%Pr=exp(phi);

%Fuel effective lower heating value
hf_0=4.3095*10^7;
deltahf=-1607.2+4.47659*T+4.00997*10^-3*T^2-6.12432*10^-7*T^3;
hf=hf_0-r*deltahf;
end
