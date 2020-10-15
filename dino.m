% simulation of the innovation/decay process
    
function exitcode=dino(logcost,rateonco,nsteps,extm,popc,flighttime,extmortcoef,nof_bins,area,disttype,ind)
    % setting the random seed    
    t1 = datetime('now','Format','dd-MMM-yyyy HH:mm:ss.SSS');
    s=second(t1,'secondofday');
    seed=s*ind; 
    seedfilename=strcat('seed.txt');
    dlmwrite(seedfilename,seed);
    rng(seed);
    
    extmort=string(extm);popchange=string(popc);
    
    % body sizes roughly follow one of the datasets shown in Lee et al. 2014
    % we track changes for this many years
    nof_bins=nof_bins+1;
    totalyears=24e7;times_in_year=linspace(0,totalyears,nof_bins);
    bodysizes_cs=getBS(nof_bins);    
        
    % steps are the distribution of mutations in cancer defences   
    mutstep_disttype=disttype; % default is string('gevrnd');
      
    % writing the parameters file
    filename_params= 'params.txt';
    params=strcat('\nnofbins: ',num2str(nof_bins),...
        '\ntotal_years: ', num2str(totalyears),...
        '\nsteps: ', mutstep_disttype,...
        '\nextrinsic mortality: ', extmort,...
        '\npopulation size: ', popchange,...
        '\nnumber of steps: ', num2str(nsteps));
    fileID = fopen(filename_params,'w');
    fprintf(fileID,params);
    fclose(fileID);
    
    cost=10.^logcost;
    
     %%%%%%
    mutrate_range_temp=logspace(-20,0,21);  
    
    % shuffling the order at which these are simulated to check the results
    % from different parts of the parameter range throughout the
    % simulations (while they are running)
    mutrate_range = mutrate_range_temp(randperm(length(mutrate_range_temp)));
    
    % get the life history traits for all the different body size scenarios
    
    % when the body size changes through macroevolutionary time
    [extrinsic_cs,LSextmort_cs,optimaldefences_cs,optimalFitness_cs,optimalLS_cs,popsize_cs]=...
    getLifeHistory(bodysizes_cs,extmort,flighttime,extmortcoef,nsteps,rateonco,...
        cost,area,popchange,times_in_year);

    % when the body size remains the same as the big ancestral size
    bodysizes_big=bodysizes_cs(1).*ones(nof_bins,1);
    [extrinsic_big,LSextmort_big,optimaldefences_big,optimalFitness_big,optimalLS_big,popsize_big]=...
    getLifeHistory(bodysizes_big,extmort,flighttime,extmortcoef,nsteps,rateonco,...
        cost,area,popchange,times_in_year);
        
    % when the body size remains the same as the small bird ancestor like
    bodysizes_small=bodysizes_cs(end).*ones(nof_bins,1);
    [extrinsic_small,LSextmort_small,optimaldefences_small,optimalFitness_small,optimalLS_small,popsize_small]=...
    getLifeHistory(bodysizes_small,extmort,flighttime,extmortcoef,nsteps,rateonco,...
        cost,area,popchange,times_in_year);
    
       
    % the file that has the last status of the simulation for all the costs
    % and all the mutation rates
    bigfilename= 'lastrows.txt';

    for j=1:length(mutrate_range)
        mutationrate=mutrate_range(j);
        for i=1:3
            if i==1
                name='changing';
                bodysizes=bodysizes_cs;                
                extrinsic=extrinsic_cs;
                LSextmort=LSextmort_cs;
                optimaldefences=optimaldefences_cs;
                optimalFitness=optimalFitness_cs;
                optimalLS=optimalLS_cs;
                popsize=popsize_cs;
            end
            if i==2
                name='big';
                bodysizes=bodysizes_big;
                extrinsic=extrinsic_big;
                LSextmort=LSextmort_big;
                optimaldefences=optimaldefences_big;
                optimalFitness=optimalFitness_big;
                optimalLS=optimalLS_big;
                popsize=popsize_big;
            end
            if i==3                
                name='small';
                bodysizes=bodysizes_small;
                extrinsic=extrinsic_small;
                LSextmort=LSextmort_small;
                optimaldefences=optimaldefences_small;
                optimalFitness=optimalFitness_small;
                optimalLS=optimalLS_small;
                popsize=popsize_small;
            end        
            % all the other dynamics happen in this function 
            % the option of having extrinsic mortality or population size
            % decoupled from the body size is only relevant for the
            % scenario where the body size is changing 
            if ~((i==2 || i==3) && (extmort~=string('changing') || popchange~=string('changing')))
                Data=dinodyn(cost, mutationrate,mutstep_disttype,nof_bins,bodysizes,...
                    times_in_year,nsteps,rateonco,num2str(j),name,...
                    extrinsic,LSextmort,optimaldefences,optimalFitness,optimalLS,popsize);
                dlmwrite(bigfilename,[i,cost,Data],'-append');
            end
        end
    end

exitcode=0;
end