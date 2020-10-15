function [extrinsic,LSextmort,optimaldefences,optimalFitness,optimalLS,popsize]=...
    getLifeHistory(bodysizes,extmort,flighttime,extmortcoef,nsteps,rateonco,...
    cost,area,population,times_in_year)
    % the function to get all the life history parameters before the
    % beginning of the macroevolutionary dynamics 
    
    % I use this function to get extrinsic mortality, see the function
    % implementation for details
    extrinsic=getExtMort(bodysizes);

    % If extrinsic mortality is decoupled from the body size, then it does not
    % change
    if extmort==string('same')
        temp=extrinsic(1);
        extrinsic=temp.*ones(length(times_in_year),1);
    % or flight can happen at a given time
    elseif extmort==string('flight')
        % take the year that flight starts
        % get the bin that flight starts
        flighthappens=find(flighttime<=times_in_year);
        flightbin=flighthappens(1);
        % change the extrinsic mortality after that time depending on the
        % extmort coef or make it the same as big animals starting from that
        % point
        extrinsic=flight(flightbin,extmortcoef,extrinsic);
    end


    % Lifespans at each timestep if extrinsic mortality was the only cause of
    % death
    % this is basically 1/extrinsic mortality, and here flight is also
    % taken into account
    LSextmort=expectedhealthylifespanExMort(extrinsic);

    % Optimal level of defences for each extrinsic mortality and each body size
    % given the cost, nsteps, and rateonco
    [optimaldefences,optimalFitness,optimalLS]=getoptimalLSdefence(extrinsic,nsteps,rateonco,bodysizes,cost);

    % the effective population size is calculated using the changes in
    % population density that is correlated with body size and an arbitrary
    % area, which here is 1e4
    popsize=getPopSize(bodysizes,area);

    % If population size is decoupled from the body size, then it does not
    % change over the course of the simulation
    if population==string('same')
        temp=popsize(1);
        popsize=temp.*ones(length(times_in_year),1);
    end
end