function meanCEL=cel(bodyS,n,mu,k,d,accuracy)
% accuracy=100 gives a pretty good approximation but feel free to use
% larger values as they're quite fast too
N = 3.72/70*10^13*bodyS;
P=linspace(1/accuracy,1-1/accuracy,accuracy); % we avoid computing at 0 and 1 but we go arbitrarily close to them when accuracy is a large number
tvalues=log(1-(1-(1-P).^(N^(-1))).^(n^(-1)))/((-1+d)*k);
meanCEL=mean(exp(-mu*tvalues));