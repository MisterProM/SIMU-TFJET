function [specific_thrust,overall_eff,C_ts,prop_eff,thermal_eff,Deltak_c,w_s,Tt0,Pt0,pi_fs,Tt2,Pt2,Tt3_prima,Pt3,pi_c,Pt4,alpha_b,Tt5_mix,Pt5,alpha_prima,T9,Pt9,pi_9,V9,V_0] = CALCULATION_PC_ID(M_0,T_0,gamma,r,alpha,P_0,e_i,Tt3_prima_val,eta_cp,e_b,Tt4,bleed,eta_tp,nozzle_velocity_coeff,e_n,M_2_val,tau_c)


Tt0=T_0*(1+(gamma-1)/2*M_0^2);
Pt0=P_0*(1+(gamma-1)/2*M_0^2)^(gamma/(gamma-1));
V_0=sqrt(2*(gamma*r/(gamma-1))*(Tt0-T_0));

pi_fs=Pt0/P_0;

%% 0-2 Intake
Tt2=T_0*(1+(gamma-1)/2*M_2_val^2);
Pt2=Pt0;

%% 2-3 Compressor

if length(M_2_val)==1
    if Tt3_prima_val==0
        Tt3_prima= tau_c*Tt2;
    else
        tau_c= Tt3_prima_val/Tt2;
    end
else
    Tt3_prima = sqrt(Tt4*T_0);
    tau_c=Tt3_prima/Tt2;
end
pi_c = tau_c^(gamma/(gamma-1));
Pt3 = Pt2*pi_c;
%% 3-4 Burner I
Pt4=Pt3*(1-e_b); 


%% 4-5 Burner II
hf_0=4.31*10^7;
if Tt3_prima_val==0
    alpha_b=((gamma*r)/(gamma-1)*(Tt4-Tt3_prima))/hf_0;
else
    alpha_b=((gamma*r)/(gamma-1)*(Tt4-Tt3_prima_val))/hf_0;
end

%% Turbine 5
if Tt3_prima_val==0
    Tt5_mix=Tt4-(Tt3_prima-Tt2);
else
    Tt5_mix=Tt4-(Tt3_prima_val-Tt2);
end
Pt5=(Tt5_mix/Tt4)^(gamma/(gamma-1))*Pt4;


%% 5-9 Nozzle
Pt9=Pt5*(1-e_n);
pi_9=Pt9/P_0;
Tt9= Tt5_mix;
T9=Tt4*(Pt3/P_0)^((1-gamma)/gamma);

V9=sqrt(2*(gamma*r/(gamma-1))*(Tt9-T9));
alpha_prima=0;

%% Performances

specific_thrust=V9-V_0;
qin = (gamma*r)/(gamma-1)*(Tt4-Tt5_mix);
Deltak_c = 0.5 * V9^2 - 0.5 * V_0^2;
w_f = specific_thrust * V_0; 
eta_th = Deltak_c / qin; 
eta_pr = w_f/Deltak_c;
if Tt3_prima_val == 0
    w_s = -(gamma*r)/(gamma-1)*((Tt4-Tt3_prima)+(Tt2-Tt0)+(T_0-T9));
else
    w_s = -(gamma*r)/(gamma-1)*((Tt4-Tt3_prima_val)+(Tt2-Tt0)+(T_0-T9));
end

prop_eff = eta_pr;
thermal_eff = eta_th;
overall_eff=eta_th*eta_pr;
C_ts=V_0/(overall_eff*hf_0);

end