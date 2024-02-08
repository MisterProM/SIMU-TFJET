function [specific_thrust,overall_eff,C_ts,prop_eff,thermal_eff,Deltak_c,W_s,Tt0,Pt0,pi_fs,Tt2,Pt2,Pt3,pi_c,Pt4,alpha_b,Tt5_mix,Pt5,alpha_prima,T9,Pt9,pi_9,V9,V_0] = CALCULATION(M_0,T_0,gamma,r,alpha,P_0,e_i,Tt3_prima_val,eta_cp,e_b,Tt4,bleed,eta_tp,nozzle_velocity_coeff,e_n)
V_0=M_0*sqrt(gamma*r*T_0);
[h,phi4_air,Cp,hf]=Functions_tables(T_0,alpha);
ht_0=h+V_0^2/2;


Tt0=T_0*(1+(gamma-1)/2*M_0^2);
[h_0,phi0,Cp_0,hf_0]=Functions_tables(Tt0,alpha);
Pt0=P_0*(1+(gamma-1)/2*M_0^2)^(gamma/(gamma-1));

pi_fs=Pt0/P_0;

%% 0-2 Intake
Tt2=Tt0;
ht_2=ht_0;
Pt2=Pt0*(1-e_i);

%% 2-3 Compressor

[ht_2,phi2,Cp_2,hf_2]=Functions_tables(Tt2,alpha);
[ht_3,phi3,Cp_3,hf_3]=Functions_tables(Tt3_prima_val,alpha);
Pt3=Pt2*(exp(phi3-phi2))^eta_cp;
pi_c = Pt3/Pt2;
%% 3-4 Burner I
Pt4=Pt3*(1-e_b); 
[ht4_air,phi4_air,Cp_4air,hf_4air]=Functions_tables(Tt4,alpha);
alpha_b=(ht4_air-ht_3)/hf_4air;

%% 4-5 Burner II
[ht_4,phi4,Cp_4,hf_4]=Functions_tables(Tt4,alpha_b);
ht5=ht_4-(ht_3-ht_2)/((1+alpha_b)*(1-bleed));

alpha_prima=alpha_b*(1-bleed);
ht5_mix=(ht5*(1+alpha_b)*(1-bleed)+ht_3*bleed)/((1+alpha_b)*(1-bleed)+bleed);
x0=0;
fun=@(x) rev_h(x(1),ht5_mix,alpha_prima);
[x, fval, exitflag] = fsolve(fun, x0);
Tt5_mix=x;
[ht_5mix,phi5mix,Cp_5mix,hf_5mix]=Functions_tables(Tt5_mix,alpha_prima);
Pt5=Pt4*(exp(phi5mix-phi4))^(1/eta_tp);
%% 5-9 Nozzle
Pt9=Pt5*(1-e_n);
pi_9=Pt9/P_0;
phi9_is=phi5mix-log(pi_9);

x1=100;
fun=@(a) rev_pr(a(1),phi9_is,alpha_prima);
[a, fval2, exitflag2] = fsolve(fun, x1);
T9_is=a;
[h_9is,phi9_is,Cp_9is,hf_9is]=Functions_tables(T9_is,alpha_prima);
h9=ht_5mix-nozzle_velocity_coeff^2*(ht_5mix-h_9is);

x3=0;
fun=@(b) rev_h(b(1),h9,alpha_prima);
[b, fval3, exitflag3] = fsolve(fun, x3);
T9=b;

V9i=sqrt(2*(ht_5mix-h_9is));
V9=V9i*nozzle_velocity_coeff;


%% Performances
specific_thrust=(1+alpha_prima)*V9-V_0;
C_ts=alpha_prima/specific_thrust; %Kg/(N*s)

qin = alpha_prima * hf;
Deltak_c = 0.5 * (1+alpha_prima) * V9^2 - 0.5 * V_0^2;
w_f = specific_thrust * V_0; 

thermal_eff = Deltak_c / qin; 
prop_eff = w_f/Deltak_c;
overall_eff=thermal_eff*prop_eff;
W_s=-((gamma*r)/(gamma-1))*((Tt4-Tt3_prima_val)+(T_0-T9));
end

