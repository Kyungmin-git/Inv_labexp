function [candidatesigma]=generate_sig2(currentsigma)

global Stepsize1 

candidatesigma=currentsigma+Stepsize1*randn;

%candidatet0=currentt0+stepsize3*randn;

if (candidatesigma<0.005) 
    candidatesigma=2*0.005-candidatesigma; end
if (candidatesigma>0.05) 
    candidatesigma=2*0.05-candidatesigma; end




end
