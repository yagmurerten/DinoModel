function L=expectedhealthylifespan(d,mu,n,bodyS,k)
% d = cancer defence
% mu = extrinsic mortality
% n = number of rate-limiting steps
% N = number of cells
% k = mutation rate per cell

% I assume a linear relationship between the number of cells and body size
% and use the human data for scaling for now (as it is the only estimate
% I could find so far)
% Moreover, here I assume all the cells in the body are equally likely to
% contribute to the oncogenic growth

% if we consider only stem cells
%N = 3.72/70*10^13*0.00002*bodyS;

N = 3.72/70*10^13*bodyS;
divisions=1000;
P=linspace(1/divisions,1-1/divisions,divisions);
tp=-log(1-(1-(1-P).^(1/N)).^(1/n))./((1-d')*k);
logL=log(sum((1-exp(-mu.*tp'))))-log(mu)-log(divisions); % log lifespan
L=exp(logL);
% We want the lifespan to be determined by extrinsic mortality when
% defences are perfect (rather than having NaN)
L(d==1)=1/mu;
