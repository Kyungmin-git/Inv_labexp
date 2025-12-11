addpath('/Users/kyuninkim/Downloads');
addpath('/Users/kyuninkim/Downloads/Pdata_inv_joint')
addpath('/Users/kyuninkim/Downloads/GA_parameter_2');
addpath('/import/c1/SOLIDEARTH/kkim39/exp_inv_store_50');
addpath('/import/c1/SOLIDEARTH/kkim39/Lib_exp');
load('Spect_store_50ml_new.mat','f_h');
addpath('/import/c1/SOLIDEARTH/kkim39/exp_inv_store_100');
addpath('/Users/kyuninkim/Downloads/Numerical_Simu');
load('spect_store_100_store2.mat','spect_store_100')

filtered_signal_raw=spect_store_100(:,1)';
filtered_signal_raw2=spect_store_100(:,2)';
global tau Stepsize1 Stepsize2 Stepsize3 Stepsize4 Stepsize5 Stepsize6
kk=22;
tau=10;
Stepsize1=5*10^(-5);
Stepsize2=2*10^(-3);
Stepsize3=1.5*10^(-7);
Stepsize4=2*10^(-4);
Stepsize5=4*10^(-7);
Stepsize6=9*10^(-4);


t0=rand(1,20)*tau;
t0_2=rand(1,20)*tau;
Q=80*1.2*10^(-6);
%qn=0.5*100*1.2*10^(-6)*ones(1,t0_num); % qn should be same size as t0
%Q=sum(qn)/tau;
sigma0=0.013;
%sigma0_2=0.015;
%qn_n=1/NN*ones(1,NN);
NN=20;
qn_var=randn(1,20)*(1/NN)*0.2;
qn_var=qn_var-mean(qn_var);
qn_d_var=randn(1,20)*(1/NN)*0.2;
qn_d_var=qn_d_var-mean(qn_d_var);
qn_n=1/NN*(ones(1,NN))+qn_var;
qn_n_d=1/NN*(ones(1,NN))+qn_d_var;


T=293.15;
R=0.02;
L_val=0.2;
L_val2=0.1;
kappa0=-8.1;
kappa0_2=-8.1;
D=0.4;
%D2=0.45;
phi=0.06;
phi2=0.45;
rho_s=2800;
Qf=20;
kappa_val=10^(kappa0);
kappa_val2=10^(kappa0_2);
m=[T,Q,R,L_val2,kappa_val,D,phi,rho_s,Qf];
m2=[T,Q,R,L_val,kappa_val2,D,phi2,rho_s,Qf];


qn=qn_n*Q*tau;
qn_d=qn_n_d*Q*tau;
qn_n_mean=mean(qn_n);

currentlpqn_n=-sum((qn_n-qn_n_mean).^2/(2*(0.05*0.2)^2));
currentlpqn_n_d=-sum((qn_n_d-qn_n_mean).^2/(2*(0.05*0.2)^2));


[Posc_in,Pspect,~,~]=forwardmodelgas15_10ps_ex2(m,t0,qn,sigma0);

[Posc_in2,Pspect2,~,~]=forwardmodelgas15_10ps_ex2(m2,t0_2,qn_d,sigma0);
Pspect(1)=0;
Pspect2(1)=0;

dt=0.0002;
t=0:dt:tau;
tlength=length(t);
currentsigma=sigma0;
currentsigma_d=sigma0;
currentkappa=kappa0;
currentgas=Q;
cur_phi=phi;
cur_gas=Q;
currentkappa2=kappa0_2;
currentqn_n=qn_n;
currentqn_n_d=qn_n_d;
cur_phi2=phi2;

currentt0=t0;
currentt0_2=t0_2;
global Icov Icov_2
Std_signal=0.5*filtered_signal_raw;
Cov_signal=diag(Std_signal);
Icov=inv(Cov_signal(10:801,10:801));
cur_gas2=Q;
Std_signal2=0.5*filtered_signal_raw2;
Cov_signal2=diag(Std_signal2);
Icov_2=inv(Cov_signal2(10:801,10:801));

  fc=5000;
    L=length(Pspect);
   
    model_amplitude_raw=abs(Pspect);
    window_size=31;
    Modelspect_raw = edgeAdjustedAverageFilter2(model_amplitude_raw, window_size);
    
    rk=filtered_signal_raw(10:801)-Modelspect_raw(10:801);
    rk_init=rk;
      fc=5000;
 
    model_amplitude_raw=abs(Pspect2);
    window_size=31;
    Modelspect_raw2 = edgeAdjustedAverageFilter2(model_amplitude_raw, window_size);
    
    rk_d=filtered_signal_raw2(10:801)-Modelspect_raw2(10:801);
rk_d_init=rk_d;
    llcurrent=(-1/2)*rk*Icov*rk'+(-1/2)*rk_d*Icov_2*rk_d';
   llcurrent_d=(-1/2)*rk_d*Icov_2*rk_d';
    llcurrent_s=(-1/2)*rk*Icov*rk';
    llcurrent_35=(-1/2)*rk*Icov*rk';
    llcurrent_50=(-1/2)*rk_d*Icov_2*rk_d';
  lpcurrent=-(cur_phi-0.006)^2/(2*0.01^2);
  lpcurrent_d=-(cur_phi2-0.9)^2/(2*0.1^2);

lpcur_add=0;
niter=2000000;
cur_Modelspect_raw=Modelspect_raw;
cur_Modelspect_raw_d=Modelspect_raw2;

kappastore=zeros(1,niter);
sigmastore=zeros(1,niter);
gasstore=zeros(1,niter);
gasstore2=zeros(1,niter);
phistore=zeros(1,niter);

kappastore2=zeros(1,niter);
phistore2=zeros(1,niter);

llstore=zeros(1,niter);
llstore_d=zeros(1,niter);
llstore_tot=zeros(1,niter);
llMAP=-inf;

t0_store=zeros(20,100);
t0_store2=zeros(20,100);

qn_n_store=zeros(20,100);
qn_n_store2=zeros(20,100);

Modelspect_raw_init=Modelspect_raw;
Modelspect_raw_d_init=Modelspect_raw2;
%}
warning('off', 'MATLAB:colon:nonIntegerIndex');
for k=1:niter
    k
    ka=rem(k,8);

    if ka==0 
    [candidatekappa,can_phi,can_gas]=generategas_3(currentkappa,cur_phi,cur_gas);
    candidatesigma=currentsigma;
    candidatet0=currentt0;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas2=cur_gas2;
    candidatesigma_d=currentsigma_d;

    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;
    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;



    elseif ka==1
    [candidatesigma]=generate_sig2(currentsigma);
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatet0=currentt0;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas2=cur_gas2;
    can_gas=cur_gas;
    candidatesigma_d=currentsigma_d;

    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;

    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;


    elseif ka==2
    [candidatet0]=generateexc_3(currentt0);
    candidatesigma=currentsigma;
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas=cur_gas;
    can_gas2=cur_gas2;
    candidatesigma_d=currentsigma_d;

    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;
    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;



     
    elseif ka==3
    [candidatekappa2,can_phi2,can_gas2]=generategas_6(currentkappa2,cur_phi2,cur_gas2);
    candidatesigma=currentsigma;
    candidatet0=currentt0;
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatet0_2=currentt0_2;
    can_gas=cur_gas;
    candidatesigma_d=currentsigma_d;
    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;
    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;


    

    elseif ka==4
    [candidatet0_2]=generateexc_3(currentt0_2);
    candidatesigma=currentsigma;
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0=currentt0;
    can_gas=cur_gas;
    can_gas2=cur_gas2;
    candidatesigma_d=currentsigma_d;

    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;

    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;


    elseif ka==5
    [candidatesigma_d]=generate_sig2(currentsigma_d);
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatesigma=currentsigma;
    candidatet0=currentt0;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas2=cur_gas2;
    can_gas=cur_gas;
    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;

    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;

    elseif ka==6
    [candidateqn_n,candidatelpqn_n]=generateqn(currentqn_n);
    candidatesigma_d=currentsigma_d;
    candidatesigma=currentsigma;
    candidatekappa=currentkappa;
    can_phi=cur_phi;
    candidatet0=currentt0;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas2=cur_gas2;
    can_gas=cur_gas;

    candidateqn_n_d=currentqn_n_d;
    candidatelpqn_n_d=currentlpqn_n_d;

    else
     [candidateqn_n_d,candidatelpqn_n_d]=generateqn(currentqn_n_d);    
    candidatekappa=currentkappa;
    candidatesigma=currentsigma;
    can_phi=cur_phi;
    candidatet0=currentt0;
    candidatekappa2=currentkappa2;
    can_phi2=cur_phi2;
    candidatet0_2=currentt0_2;
    can_gas2=cur_gas2;
    can_gas=cur_gas;
    candidateqn_n=currentqn_n;
    candidatelpqn_n=currentlpqn_n;
    candidatesigma_d=currentsigma_d;



    end
    %sigma=candidatesigma;
    %kappa=candidatekappa;
    %t0=candidatet0;
    
    kappa_val=10^(candidatekappa);
    kappa_val2=10^(candidatekappa2);

    candidateqn=candidateqn_n*can_gas*tau;
    candidateqn_d=candidateqn_n_d*can_gas2*tau;

    m=[T,can_gas,R,L_val2,kappa_val,D,can_phi,rho_s,Qf];
    m2=[T,can_gas2,R,L_val,kappa_val2,D,can_phi2,rho_s,Qf];
    if ka==0 || ka==1 || ka==2 || ka==6
    


   
    [Posc_in,Pspect,~,~]=forwardmodelgas15_10ps_ex2(m,candidatet0,candidateqn,candidatesigma);
    Pspect(1)=0;
    %[Posc_in2,Pspect2,~,~]=forwardmodelgas15_10ps_ex2(m2,candidatet0_2,qn,candidatesigma_d);
    
     
    %Posc_spect=abs(Posc_can);
    fc=5000;
  
    model_amplitude_raw=abs(Pspect);
    %model_amplitude_raw_d=abs(Pspect2);
    window_size=31;
    Modelspect_raw = edgeAdjustedAverageFilter2(model_amplitude_raw, window_size);
    %Modelspect_raw_d=edgeAdjustedAverageFilter2(model_amplitude_raw_d,window_size);
    rk=(filtered_signal_raw(10:801)-Modelspect_raw(10:801));
    rk2=(filtered_signal_raw2(10:801)-Modelspect_raw2(10:801));
        fc=5000;
    L=length(Posc_in2);
    df1=(1/dt)/length(Posc_in2);
   
    L2=round(L/2+1);
    llcandidate_s=(-1/2)*rk*Icov*rk';
    llcandidate_d=llcurrent_d;
    llcandidate=llcandidate_s+llcandidate_d;

    sigma_diff=candidatesigma-candidatesigma_d;
   
    lpcandidate=-(can_phi-0.006)^2/(2*0.01^2)-(sigma_diff)^2/(2*0.0005^2);
    lpcandidate_d=-(can_phi2-0.9)^2/(2*0.1^2);

    logalpha= llcandidate+lpcandidate+lpcandidate_d+candidatelpqn_n+candidatelpqn_n_d...
        -llcurrent-lpcurrent-lpcurrent_d-currentlpqn_n-currentlpqn_n_d;

    
    else
        %[Posc_in,Pspect,~,~]=forwardmodelgas15_10ps_ex2(m,candidatet0,candidateqn,candidatesigma);
        [Posc_in2,Pspect2,~,~]=forwardmodelgas15_10ps_ex2(m2,candidatet0_2,candidateqn_d,candidatesigma_d);
   Pspect2(1)=0;
%model_amplitude_raw=abs(Pspect);
    model_amplitude_raw_d=abs(Pspect2);
    window_size=31;
    %Modelspect_raw=edgeAdjustedAverageFilter2(model_amplitude_raw,window_size);
    Modelspect_raw2 = edgeAdjustedAverageFilter2(model_amplitude_raw_d, window_size);
    rk=(filtered_signal_raw(10:801)-Modelspect_raw(10:801));
    rk2=(filtered_signal_raw2(10:801)-Modelspect_raw2(10:801));
    
    llcandidate_d=(-1/2)*rk2*Icov_2*rk2';
    llcandidate_s=llcurrent_s;
    llcandidate=llcandidate_d+llcandidate_s;
      sigma_diff=candidatesigma-candidatesigma_d;
   
    lpcandidate=-(can_phi-0.006)^2/(2*0.01^2)-(sigma_diff)^2/(2*0.0005^2);
    lpcandidate_d=-(can_phi2-0.9)^2/(2*0.1^2);

    logalpha= llcandidate+lpcandidate+lpcandidate_d+candidatelpqn_n+candidatelpqn_n_d...
        -llcurrent-lpcurrent-lpcurrent_d-currentlpqn_n-currentlpqn_n_d;


    end

     if (logalpha>0)
         logalpha=0;
     end
    
    logt=log(rand());

    if (logt<logalpha)
       currentt0=candidatet0;
       currentt0_2=candidatet0_2;
       currentkappa=candidatekappa;
       currentsigma=candidatesigma;
       currentsigma_d=candidatesigma_d;
       %currentgas=candidategas;
       cur_phi=can_phi;
       currentkappa2=candidatekappa2;
       cur_phi2=can_phi2;
       cur_gas=can_gas;
       cur_gas2=can_gas2;
       lpcurrent=lpcandidate;
       llcurrent=llcandidate;
       llcurrent_s=llcandidate_s;
       llcurrent_d=llcandidate_d;
       currentqn_n=candidateqn_n;
       currentlpqn_n=candidatelpqn_n;
       currentqn_n_d=candidateqn_n_d;
       currentlpqn_n_d=candidatelpqn_n_d;
       

       acceptance(k)=1;
       cur_Modelspect_raw=Modelspect_raw;
       cur_Modelspect_raw_d=Modelspect_raw2;

    else
       
       acceptance(k)=0;
    end
       
       
     %value store
     kappastore(1,k)=currentkappa;
     sigmastore(1,k)=currentsigma;
     sigmastore2(1,k)=currentsigma_d;
     phistore(1,k)=cur_phi;
     gasstore(1,k)=cur_gas;
     gasstore2(1,k)=cur_gas2;
     kappastore2(1,k)=currentkappa2;
     phistore2(1,k)=cur_phi2;
     %gasstore(1,k)=currentgas;

     llstore(1,k)=llcurrent;

     %store MAP solution
     if llcandidate>llMAP
     llMAP=llcandidate;
     kappa_MAP=currentkappa;
     sigma_MAP=currentsigma;
     sigma_MAP2=currentsigma_d;
     kappa_MAP2=currentkappa2;
    
     phi_MAP=cur_phi;
     phi_MAP2=cur_phi2;
     %gas_MAP=currentgas;
     t0_MAP=currentt0;
     t0_MAP2=currentt0_2;
     qn_n_MAP=currentqn_n;
     qn_n_MAP2=currentqn_n_d;
     spect_MAP=Modelspect_raw;
     spect_MAP2=Modelspect_raw2;
     end

     if mod(k,20000)==0 || k==100 || k==1000


         h13=figure(13);
         [~,edges] =histcounts(log10(10.^(kappastore(1:k))),50);
         histogram(10.^(kappastore(1:k)),10.^edges,'Normalization','Probability');
         xlim([10^(-10) 10^(-6)])
         xticks([10^(-9) 10^(-8) 10^(-7)])
         xline(10^kappa_MAP,'r')
         xlabel('m^{2}')
         title('Permeability')
         set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h13,sprintf('MCMC_source_%d.fig',kk),'compact')
         saveas(h13,sprintf('MCMC_source_%d.jpg',kk))
      h14=figure(14);
         
         histogram(sigmastore(1:k),'Normalization','Probability','NumBins',50);
         xlim([0 0.05])
         xticks([0 0.02 0.04])
         xline(sigma_MAP,'r')
         xlabel('time (s)')
         title('Sigma')
         %set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h14,sprintf('MCMC_source2_%d.fig',kk),'compact')
         saveas(h14,sprintf('MCMC_source2_%d.jpg',kk))
     
             
         h17=figure(17);
         
         histogram(100*phistore(1:k),'Normalization','Probability','NumBins',50);
         xlim([0 100])
         xline(100*phi_MAP,'r')
         %xticks([0 0.02 0.04])
         xlabel('porosity (%)')
         title('Sigma')
         %set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h17,sprintf('MCMC_source4_%d.fig',kk),'compact')
         saveas(h17,sprintf('MCMC_source4_%d.jpg',kk))

     lf=length(f_h);
     h15=figure(15);
     plot(f_h,filtered_signal_raw)
     hold on
     plot(f_h,spect_MAP(1:lf))
     hold off
     set(gca,'FontSize',18)
     xlim([0 60])
     xlabel('Frequency (Hz)')
     ylabel('Pressure (Pa*s)')
     saveas(h15,sprintf('MCMC_source5_%d.jpg',kk))

     
     h18=figure(18);
         [~,edges] =histcounts(log10(10.^(kappastore2(1:k))),50);
         histogram(10.^(kappastore2(1:k)),10.^edges,'Normalization','Probability');
         xlim([10^(-10) 10^(-6)])
         xticks([10^(-9) 10^(-8) 10^(-7)])
         xline(10^kappa_MAP2,'r')
         xlabel('m^{2}')
         title('Permeability')
         set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h18,sprintf('MCMC_source_50_%d.fig',kk),'compact')
         saveas(h18,sprintf('MCMC_source_50_%d.jpg',kk))

         h21=figure(21);
         
         histogram(gasstore(1:k)/(1.2*10^(-6)),'Normalization','Probability','NumBins',50);
         xlim([0 200])
        
         %xticks([0 0.02 0.04])
         xlabel('gas flow rate (ml/s)-9 hole')
         title('')
         %set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h21,sprintf('MCMC_source6_%d.fig',kk),'compact')
         saveas(h21,sprintf('MCMC_source6_%d.jpg',kk))

                h22=figure(22);
         
         histogram(100*phistore2(1:k),'Normalization','Probability','NumBins',50);
         xlim([0 100])
         xline(100*phi_MAP2,'r')
         %xticks([0 0.02 0.04])
         xlabel('porosity (%)')
         title('Sigma')
         %set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h22,sprintf('MCMC_source7_%d.fig',kk),'compact')
         saveas(h22,sprintf('MCMC_source7_%d.jpg',kk))


        h19=figure(19);
         
         histogram(gasstore2(1:k)/(1.2*10^(-6)),'Normalization','Probability','NumBins',50);
         xlim([0 200])
         
         %xticks([0 0.02 0.04])
         xlabel('gas flow rate - dense')
         title('')
         %set(gca,'XScale','log')
         set(gca,'YTickLabel',[])
         savefig(h19,sprintf('MCMC_source8_%d.fig',kk),'compact')
         saveas(h19,sprintf('MCMC_source8_%d.jpg',kk))


             h20=figure(20);
     plot(f_h,filtered_signal_raw2)
     hold on
     plot(f_h,spect_MAP2(1:lf))
     hold off
     set(gca,'FontSize',18)
     xlim([0 60])
     xlabel('Frequency (Hz)')
     ylabel('Pressure (Pa*s)')
     saveas(h20,sprintf('MCMC_fit8_%d.jpg',kk))

     accept1=acceptance(1:8:k);
     accept2=acceptance(2:8:k);
     accept3=acceptance(3:8:k);
     accept4=acceptance(4:8:k);
     accept5=acceptance(5:8:k);
     accept6=acceptance(6:8:k);
     accept7=acceptance(7:8:k);
     accept10=acceptance(8:8:k);

     acceptance1=sum(accept1)/length(accept1);
     acceptance2=sum(accept2)/length(accept2);
     acceptance3=sum(accept3)/length(accept3);
     acceptance4=sum(accept4)/length(accept4);
     acceptance5=sum(accept5)/length(accept5);
     acceptance6=sum(accept6)/length(accept6);
     acceptance7=sum(accept7)/length(accept7);
     acceptance10=sum(accept10)/length(accept10);


     save(sprintf('expinvresult_%d.mat',kk),'-v7.3')

      




     end

     if mod(k,20000)==0 
    
     t0_store(:,k/20000)=currentt0';
     t0_store2(:,k/20000)=currentt0_2';

     qn_n_store(:,k/20000)=currentqn_n';
     qn_n_store2(:,k/20000)=currentqn_n_d';
     else
     end
   






end
