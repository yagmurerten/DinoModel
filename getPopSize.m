function m=getPopSize(bodS,area)

% function to find the effective population size
% by plugging in body size and area (to later change this to see what 
% effect it has)

% using population density and body size correlation for birds from 
% Juanes 1986
popDensity=10.^(-0.49.*log10(1000.*bodS)+1.96);

popSize=area.*popDensity;
m=popSize';