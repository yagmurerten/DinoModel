function m=getExtMort(bodS)

% we get the extrinsic mortality scaling from McCarthy et al. 2008 
% termed instantaneous mortality rate there
% there they report a mean of b=-0.21 in their posterior distribution
% and we get ln(a)=-1.8 using plotdigitizer from their appendices

%therefore we have
m=exp(-1.8).*bodS.^(-0.21);
