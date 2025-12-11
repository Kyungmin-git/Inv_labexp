function [candidatet0]=generateexc_3(currentt0)
global tau
%Stepsize=0.02*abs(normrnd(0,1));

Ixshift=randperm(20,1);
currentt0(Ixshift)=currentt0(Ixshift)+normrnd(0,0.03,1);
currentt0=abs(currentt0);
idxt0=find(currentt0>tau);
currentt0(idxt0)=2*tau-currentt0(idxt0);
candidatet0=currentt0;



end
