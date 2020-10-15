function [optd,fitness,lifespan]=getoptimalLSdefence(mu,n,k,bodysize,cost)
    accuracy=1e2;
    % very high accuracy crashes matlab
    if accuracy > 1e5
        accuracy = 1e5;
    end
    d=linspace(0,1,accuracy);
    optd=zeros(length(bodysize),1);fitness=zeros(length(bodysize),1);lifespan=zeros(length(bodysize),1);
    for i=1:length(bodysize)
        lifespans=expectedhealthylifespan(d,mu(i),n,bodysize(i),k);
        if cost~=0
            costcoef=1-d.^(1/cost); 
            % we calculate fitness taking the trade-off with the defences
            % into account. The strength of this trade-off depends on the
            % cost
            scaledLS=costcoef.*lifespans; 
            [maxfit,index]=max(scaledLS);
            ex_d=d(index);ex_fit=maxfit;ex_LS=lifespans(index);
            % to get a better estimate of the optimal fitness, I zoom in at
            % this interval
            dtemp=linspace(d(index)-1/accuracy,d(index)+1/accuracy,1e3);
            dtemp(dtemp>1.0)=1.0;dtemp(dtemp<0)=0.0;
            lifespans_temp=expectedhealthylifespan(dtemp,mu(i),n,bodysize(i),k);
            costcoef=1-dtemp.^(1/cost);
            % scaledLS is essentially 'fitness', scaled by the cost
            % coefficient
            scaledLS_temp=costcoef.*lifespans_temp;
            [maxfit,index]=max(scaledLS_temp);
            optd(i)=dtemp(index);fitness(i)=maxfit;lifespan(i)=lifespans_temp(index);
        else
            % no cost -> perfect defences, and the lifespan=fitness=inverse
            % of extrinsic mortality rate
            optd(i)=1.0;fitness(i)=1/mu(i);lifespan(i)=1/mu(i);
        end
        
    end
end