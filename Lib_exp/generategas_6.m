function [candidatekappa,can_phi,can_gas]=generategas_6(currentkappa,cur_phi,cur_gas)

global Stepsize2 Stepsize6 Stepsize5

%candidatesigma=currentsigma+Stepsize1*randn;
candidatekappa=currentkappa+Stepsize2*randn;

can_phi=cur_phi+Stepsize6*randn;

can_gas=cur_gas+Stepsize5*randn;
%candidatet0=currentt0+stepsize3*randn;


if (candidatekappa<-11) 
    candidatekappa=2*(-11)-candidatekappa; end
if (candidatekappa>-6.3) 
    candidatekappa=2*(-6.3)-candidatekappa; end

if (can_phi<10^(-3)) 
    can_phi=2*(10^(-3))-can_phi; end
if (can_phi>(1)) 
    can_phi=2*(1)-can_phi; end

if (can_gas<0) 
    can_gas=-can_gas; end
if (can_gas>(6*10^(-4))) 
    can_gas=2*(6*10^(-4))-can_gas; end
