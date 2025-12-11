function [candidateqn_n,candidatelpqn_n]=generateqn(currentqn_n)

%Stepsize=0.02*abs(normrnd(0,1));
NN=20;
Ixshift=randperm(20,1);
currentqn_n(Ixshift)=currentqn_n(Ixshift)+0.05*normrnd(0,0.08,1);
currentqn_n=abs(currentqn_n);
candidateqn_n=currentqn_n/sum(currentqn_n);
qn_n_mean=1/NN*(ones(1,NN));
candidatelpqn_n=-sum((candidateqn_n-qn_n_mean).^2/(2*(0.05*0.2)^2));




end
