function [specific_thrust,overall_eff,C_ts,prop_eff,thermal_eff,Deltak_c,w_s,Tt0,Pt0,pi_fs,Tt2,Pt2,Tt25,Pt25,Tt25_p,Pt25_p,Tt3_prima,Pt3_prima,pi_c,Pt4,alpha_b,Tt5_mix,Pt5,alpha_prima,T9,Pt9,pi_9,V9,V_0,s_0,s0,s2,s25,s25_p,s3,s4,s5,s9] = CALCULATION_IC_TT3P(M_0,T_0,gamma,r,alpha,P_0,e_i,Tt3_val_prima,Delta_T_Pc_values,eta_cp,e_b,Tt4,bleed,eta_tp,nozzle_velocity_coeff,e_n,e_interp,diff_T,Tt25_val)
V_0=M_0*sqrt(gamma*r*T_0);
[h,phi4_air,Cp,hf]=Functions_tables(T_0,alpha);
ht_0=h+V_0^2/2;
s_0 = r*phi4_air-r*log(P_0/(exp(phi4_air)));

Tt0=T_0*(1+(gamma-1)/2*M_0^2);
Pt0=P_0*(1+(gamma-1)/2*M_0^2)^(gamma/(gamma-1));
[h_0,phi0,Cp_0,hf_0]=Functions_tables(Tt0,alpha);
s0 = r*phi0-r*log(Pt0/(exp(phi0)));

pi_fs=Pt0/P_0;

%% 0-2 Intake
Tt2=Tt0;
Pt2=Pt0*(1-e_i);

%% 2-3 Compressor
% TO MAXIMIZE WORK WITH INTERCOOLING
if diff_T == 0
    Tt25_p = Tt0;
    
    Tt25 = nthroot(Tt25_p,3)*nthroot(Tt4,3)*nthroot(T_0,3);
    
    Tt3_prima = nthroot(Tt25_p,3)*nthroot(Tt4,3)*nthroot(T_0,3);
    
    [ht_2,phi2,Cp_2,hf_2]=Functions_tables(Tt2,alpha);
    s2 = r*phi2-r*log(Pt2/(exp(phi2)));
    
    if Tt25_val==0
        [ht_25,phi25,Cp_25,hf_25]=Functions_tables(Tt25,alpha);
    else
        [ht_25,phi25,Cp_25,hf_25]=Functions_tables(Tt25_val,alpha);
    end
    
    Pt25=Pt2*(exp(phi25-phi2))^sqrt(eta_cp);
    s25 = r*phi25-r*log(Pt25/(exp(phi25)));
    
    Pt25_p=Pt25*(1-e_interp); %CHECK IF IT IS ASSOCIATED TO THE TEMPERATURE THAT WANTS TO BE COOLED
    [ht_25_p,phi25_p,Cp_25_p,hf_25_p]=Functions_tables(Tt25_p,alpha);
    s25_p = r*phi25_p-r*log(Pt25_p/(exp(phi25_p)));
    
    if Tt3_val_prima == 0
        [ht_3p,phi3,Cp_3,hf_3]=Functions_tables(Tt3_prima,alpha);
    else
        [ht_3p,phi3,Cp_3,hf_3]=Functions_tables(Tt3_val_prima,alpha);
    end
else
    [ht_2,phi2,Cp_2,hf_2]=Functions_tables(Tt2,alpha);
    [ht_25,phi25,Cp_25,hf_25]=Functions_tables(Tt25_val,alpha);
    Pt25=Pt2*(exp(phi25-phi2))^sqrt(eta_cp);
    Pt25_p=Pt25*(1-e_interp);
    Tt25_p = (Tt25_val-Tt0)*diff_T + Tt0;
    [ht_25_p,phi25_p,Cp_25_p,hf_25_p]=Functions_tables(Tt25_p,alpha);
    [ht_3p,phi3,Cp_3,hf_3]=Functions_tables(Tt3_val_prima+(Tt25_val-Tt0)*diff_T,alpha);
end
Pt3_prima=Pt25_p*(exp(phi3-phi25_p))^sqrt(eta_cp);


s3 = r*phi3-r*log(Pt3_prima/(exp(phi3)));
pi_c_1=Pt25/Pt2;
pi_c_2 = Pt3_prima/Pt25_p;
pi_c=pi_c_1*pi_c_2;
%% 3-4 Burner I
Pt4=Pt3_prima*(1-e_b); 
[ht4_air,phi4_air,Cp_4air,hf_4air]=Functions_tables(Tt4,alpha);

alpha_b=(ht4_air-ht_3p)/hf_4air;

%% 4-5 Burner II
[ht_4,phi4,Cp_4,hf_4]=Functions_tables(Tt4,alpha_b);
s4 = r*phi4-r*log(Pt4/(exp(phi4)));
ht5=ht_4-((ht_3p-ht_25_p)+(ht_25-ht_2))/((1+alpha_b)*(1-bleed));

alpha_prima=alpha_b*(1-bleed);
ht5_mix=(ht5*(1+alpha_b)*(1-bleed)+ht_3p*bleed)/((1+alpha_b)*(1-bleed)+bleed);
x0=0;
fun=@(x) rev_h(x(1),ht5_mix,alpha_prima);
[x, fval, exitflag] = fsolve(fun, x0);
Tt5_mix=x;
[ht_5mix,phi5mix,Cp_5mix,hf_5mix]=Functions_tables(Tt5_mix,alpha_prima);
Pt5=Pt4*(exp(phi5mix-phi4))^(1/eta_tp);
s5 = r*phi5mix-r*log(Pt5/(exp(phi5mix)));
%% 5-9 Nozzle
Pt9=Pt5*(1-e_n);
pi_9=Pt9/P_0;
phi9_is=phi5mix-log(pi_9);

x1=10;
fun=@(a) rev_pr(a(1),phi9_is,alpha_prima);
[a, fval2, exitflag2] = fsolve(fun, x1);
T9_is=a;
[h_9is,phi9_is,Cp_9is,hf_9is]=Functions_tables(T9_is,alpha_prima);
h9=ht_5mix-nozzle_velocity_coeff^2*(ht_5mix-h_9is);
x3=0;
fun=@(b) rev_h(b(1),h9,alpha_prima);
[b, fval3, exitflag3] = fsolve(fun, x3);
T9=b;

[h_9,phi9,Cp_9,hf_9]=Functions_tables(T9,alpha_prima);
s9 = r*phi9-r*log(Pt9/(exp(phi9)));

ht9 = ht_5mix+eta_tp*(h_9is-ht_5mix);

V9i=sqrt(2*(ht_5mix-h_9is));
V9=V9i*nozzle_velocity_coeff;








%% Performances

    specific_thrust=(1+alpha_prima)*V9-V_0;
    qin = ht_4-ht_3p;
    Deltak_c = 0.5 * (1+alpha_prima) * V9^2 - 0.5 * V_0^2; 
    eta_th = Deltak_c / qin; 
    eta_pr = specific_thrust * V_0/Deltak_c;
    
    C_ts=alpha_prima/specific_thrust; %g/kN*s
    
    prop_eff = eta_pr;
    thermal_eff = eta_th;
    overall_eff=eta_th*eta_pr;
    w_s=-((ht_4-ht_3p)+(ht_25_p-ht_25)+(h-h_9));
end