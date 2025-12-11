function [Posc,Pspect,A_res,A_exc]=forwardmodelgas15_10ps_ex2(m,t0,qn,sigma_k)

%% PHYSICAL PARAMETERS

mu_g=1e-5;                                                                  %gas viscosity
M=0.029;                                                                    %molecular weight of gas (water vapor)
Rg=8.3145;  
Pex=101325;   

T=m(1);
Q=m(2);
R=m(3);
L=m(4);
kappa=m(5);
D=m(6);
phi=m(7);

S=pi*R^2;
%% AUXILIARY PARAMETERS
beta_a=S*phi*M/(Rg*T);
beta_b=mu_g*phi/(kappa*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_c=Pex*M/(Rg*T*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_d=2*Q/(S*phi*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));
beta_e=S*M*D/(Rg*T);
P0=Pex+mu_g*Rg*T*Q*L/(S*kappa*M*(Pex-Rg*T*Q^2/(S^2*phi^2*M*Pex)));

%% COEFFICIENTS OF THE HARMONIC OSCILLATOR
GAMMA0=1;
GAMMA1=(2*(beta_a *beta_d+beta_b *beta_e )*L+beta_a *beta_b *L^2)/(2*beta_a );
GAMMA2=(2*beta_c *beta_e *L+beta_a *beta_c *L^2)/(2*beta_a );
gamma0=beta_b*L/beta_a;
gamma1=beta_c*L/beta_a;

%% NATURAL FREQUENCY AND CRITICAL THICKNESS
fn=sqrt((sqrt((GAMMA2*gamma0^2+gamma1^2)^2-GAMMA1^2*gamma0^2*gamma1^2)-GAMMA2*gamma0^2)/(GAMMA2*gamma1^2))/(2*pi);
a0=(beta_a*beta_b*L^2+2*beta_a*beta_d*L)^2*beta_b^2-4*beta_a^2*beta_b^2*beta_c*L^2-4*beta_a^2*beta_c^2;
a1=4*(beta_a*beta_b*L^2+2*beta_a*beta_d*L)*beta_b^3*L-8*beta_a*beta_b^2*beta_c*L;
a2=4*beta_b^4*L^2;
Dcrit=(Rg*T/(S*M))*(-a1+sqrt(a1^2-4*a0*a2))/(2*a2);
%disp(['Natural frequency of the oscillator is: ' num2str(fn) ' Hz'])
%disp(['Critical thickness of the gas pocket: ' num2str(Dcrit) ' m'])
%% Simulation time unit
global tau

max_freq=5000;
dt=1/max_freq;
%NN=LPendt-LPstartt;
%inter=1/tau;
%t=0:dt:tau; 
time1=0:dt:tau;
df=(1/dt)/length(time1);
Lw=length(time1);
%% SUPPLY OF VOLATILES TO THE GAS POCKET  

    %qn_ave=Q*tau/length(t0);
    %qn=qn_ave.*qn_n;
%% CALCULATION OF THE PRESSURE EVOLUTION IN THE GAS POCKET AND GROUND DISPLACEMENT IN THE FREQUENCY DOMAIN [using equation (C11)]
% note that from the frequency domain we calculate the steady-state pressure
% evolution, and thus we do not account for the transient state. That
% is why the pressure time series calculated with this approach does not
% exactly coincide with the pressure time series calculated directly in the time domain 
% (equation (15)) during the first seconds of simulation. 
mm=1;
i=1i;
global distance
u_z=zeros(1,Lw);
A_p=zeros(1,Lw);
A_res=zeros(1,Lw);
A_exc=zeros(1,Lw);
A_path=zeros(1,Lw);
N_exc=length(qn);


i=1i;
for j=0:1/tau:1/(2*dt) % calculate only positive frequency
    ome=2*pi*j;
    freq_teo(mm)=j; 

    %source
    A_res(mm)=(gamma0*(i*ome)^0+gamma1*(i*ome)^1)/(GAMMA0*(i*ome)^0+GAMMA1*(i*ome)^1+GAMMA2*(i*ome)^2);  
    A_exc(mm)=sum(qn.*exp(-i*ome*t0).*exp(-(1/2)*sigma_k^2.*ome^2));                                               
    A_p(mm)=A_res(mm)*A_exc(mm);                                                    


 

    mm=mm+1;
end

% reconstruction of the signals in the time domain
P=ifft(A_p,'symmetric')*sqrt(Lw)/sqrt(dt/df);
Pspect=A_p;
Posc=P+Pex-P0;




    %disp(['Natural frequency of the oscillator is: ' num2str(fn) ' Hz'])
    %disp(['Critical thickness of the gas pocket: ' num2str(Dcrit) ' m'])
%% Calculate Synthetic amplitude
