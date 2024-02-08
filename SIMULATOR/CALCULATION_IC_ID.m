function [specific_thrust,overall_eff,C_ts,prop_eff,thermal_eff,Deltak_c,w_s,Tt0,Pt0,pi_fs,Tt2,Pt2,Tt25,Pt25,Tt25_p,Pt25_p,Tt3_prima,Pt3,pi_c,Pt4,alpha_b,Tt5_mix,Pt5,alpha_prima,T9,Pt9,pi_9,V9,V_0] = CALCULATION_IC_ID(M_0,T_0,gamma,r,alpha,P_0,e_i,Tt3_prima_val,Delta_T_Pc_values,eta_cp,e_b,Tt4,bleed,eta_tp,nozzle_velocity_coeff,e_n,e_interp,diff_T,Tt25_p_val)


Tt0=T_0*(1+(gamma-1)/2*M_0^2);
Pt0=P_0*(1+(gamma-1)/2*M_0^2)^(gamma/(gamma-1));
V_0=sqrt(2*(gamma*r/(gamma-1))*(Tt0-T_0));

pi_fs=Pt0/P_0;

%% 0-2 Intake
Tt2=Tt0;
Pt2=Pt0;

%% 2-3 Compressor
% TO MAXIMIZE WORK WITH INTERCOOLING
if Tt25_p_val==0
    Tt25_p = Tt0;
else
    Tt25_p =Tt25_p_val;
end

Tt25 = nthroot(Tt25_p,3)*nthroot(Tt4,3)*nthroot(T_0,3);
Tt3_prima = nthroot(Tt25_p,3)*nthroot(Tt4,3)*nthroot(T_0,3);

Pt25=(Tt0/Tt25)^(gamma/(1-gamma))*Pt0;
Pt25_p=Pt25;
if Tt3_prima_val==0
    Pt3=(Tt25_p/Tt3_prima)^(gamma/(1-gamma))*Pt25_p;
else
    Pt3=(Tt25_p/Tt3_prima_val)^(gamma/(1-gamma))*Pt25_p;
end

pi_c=(Pt3/Pt25_p)*(Pt25/Pt2);
%% 3-4 Burner I
Pt4=Pt3; 


%% 4-5 Burner II
hf_0=4.31*10^7;
if Tt3_prima_val==0
    alpha_b=((gamma*r)/(gamma-1)*(Tt4-Tt3_prima))/hf_0;
else
    alpha_b=((gamma*r)/(gamma-1)*(Tt4-Tt3_prima_val))/hf_0;
end

%% Turbine 5
if Tt3_prima_val==0
    Tt5_mix=Tt4-((Tt3_prima-Tt25_p)+(Tt25-Tt2));
else
    Tt5_mix=Tt4-((Tt3_prima_val-Tt25_p)+(Tt25-Tt2));
end

Pt5=(Tt5_mix/Tt4)^(gamma/(gamma-1))*Pt4;


%% 5-9 Nozzle
Pt9=Pt5;
pi_9=Pt9/P_0;
Tt9= Tt5_mix;
T9=Tt4*(Pt3/P_0)^((1-gamma)/gamma);

V9=sqrt(2*(gamma*r/(gamma-1))*(Tt9-T9));
alpha_prima=0;

%% Performances

specific_thrust=V9-V_0;
qin = (gamma*r)/(gamma-1)*(Tt4-Tt3_prima);
Deltak_c = 0.5 * V9^2 - 0.5 * V_0^2;
w_f = specific_thrust * V_0; 
eta_th = Deltak_c / qin; 
eta_pr = w_f/Deltak_c;
if Tt3_prima_val == 0
    w_s = -(gamma*r)/(gamma-1)*(Tt4-Tt3_prima+T_0-Tt4*(Tt25_p*T_0)/(Tt3_prima*Tt25)+Tt25_p-Tt25);
else
    w_s = -(gamma*r)/(gamma-1)*(Tt4-Tt3_prima_val+T_0-Tt4*(Tt25_p*T_0)/(Tt3_prima_val*Tt25)+Tt25_p-Tt25);
end

prop_eff = eta_pr;
thermal_eff = eta_th;
overall_eff=eta_th*eta_pr;
C_ts=V_0/(overall_eff*hf_0);

end