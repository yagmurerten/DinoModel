function DataEnd=dinodyn(cost,mutationrate,mutstep_disttype,nof_bins,bodysizes,times_in_year,...
    nsteps,rateonco,fileindex,name,extrinsic,LSextmort,optimaldefences,optimalFitness,optimalLS,popsize)

filename= strcat('dinodata_',name,'_',fileindex,'.txt');
mutfixfile=strcat('distmuts_',name,'_',fileindex,'.txt');

% initalization
currentbin=1;
mutfixes=[0 0]; % this tracks how many positive and negative mutations have fixed
costnew=1; costold=1; 
% cost multiplier is equal to 1 if there are no costs
% and depends on the defence level if there are costs (calculated later)

currenttime=0.;


% initial defence, we assume they start with the optimal defences
defence=optimaldefences(1);
accuracy=10000; d=linspace(0,1,accuracy);

if mutstep_disttype==string('gevrnd')
    K=-0.75; sigma=1; mu=-1;
    mutations=max(-1,gevinv(0.0005:0.001:0.9995,K,sigma,mu)); % this gives 100 mutations, each of them is equally likely; 
elseif mutstep_disttype==string('gaussian')
    mean=0.;sd=0.3;
    mutations=max(-1,norminv(0.0005:0.001:0.9995,mean,sd)); 
end

newbin=1; % this tells the code to recalculate lifespans, etc., as we're in a new bin
%%%% new function
while currenttime<times_in_year(end)
    if newbin
        currentpopsize=popsize(currentbin);
        currentbodysize=bodysizes(currentbin);
        currentextmort=extrinsic(currentbin);
        current_od=optimaldefences(currentbin);
        current_oFitness=optimalFitness(currentbin);
        current_oLifespan=optimalLS(currentbin);
        current_LS_extmort=LSextmort(currentbin);
        LScancer=expectedhealthylifespanCancer(defence,nsteps,currentbodysize,rateonco);
        cancerp=cel(currentbodysize,nsteps,currentextmort,rateonco,defence,1000);
        lifespan_d=expectedhealthylifespan(d,currentextmort,nsteps,currentbodysize,rateonco); % lifespans for all d values
        currentlifespan=interp1(d,lifespan_d, defence);
        if cost==0
            costcoef=1; % if there is no cost, no multiplier for the lifespans in the next eq.
        elseif cost~=0
            costcoef=(1-defence^(1/cost)); % mutiplier with the defence
        end
        % current fitness depends on the lifespan and how costly are the
        % defences
        currentfitness=costcoef*currentlifespan;
        Data=[currenttime defence currentbodysize currentpopsize currentextmort currentlifespan currentfitness mutationrate nsteps rateonco current_LS_extmort LScancer current_od current_oFitness current_oLifespan cancerp mutfixes(1)]; 
        dlmwrite(filename,Data,'-append'); % I write data in each new bin to lower the memory use
        newbin=0; % job done
    end
    
    % getting the changes in defences for all possible different mutations
    mutated_defence=(mutations<0).*(defence-defence*(-mutations))+(mutations>=0).*(defence+mutations*(1-defence));
    lifespan_current=interp1(d,lifespan_d,defence);
    lifespan_mutated=interp1(d,lifespan_d,mutated_defence);
    costmultiplier_current=(1-defence^(1/cost));
    costmultiplier_mutated=(1-mutated_defence.^(1/cost));
    selection=((lifespan_mutated.*costmultiplier_mutated)-(lifespan_current*costmultiplier_current))/(lifespan_current*costmultiplier_current);
    fixed=0;
	rate_s = (4*mutationrate/(length(mutations))*currentpopsize.*selection)./(1-exp(-4*currentpopsize.*selection)); % NOTE THIS IS ASSUMING RPECISELY 100 CATEGORIES IN GEVINV
    rate_s(selection==0)=mutationrate;
    % these mutations 'compete'
    [gen_till_event,which_mutation_fixed]=min(exprnd(1./rate_s));% how many generations until the next fixation; this chooses the shortest fixation
    t_till_event=gen_till_event*lifespan_current; % time until this event (using lifespan that was valid. Could also use mean of old and new lifespan to approximate this...)
    try
            bin=find(currenttime+t_till_event>times_in_year); bin=bin(end);
    catch
            currenttime;
            lifespans(2);
    end
%     bin=find(currenttime+t_till_event>times_in_year); bin=bin(end);
    if bin>nof_bins
        disp('End of simu')
        break
    else
        if bin==currentbin % we didn't move to a new bin, so the 'plus' event really did happen first
            currenttime=currenttime+t_till_event;
            fixed=1;
            defence=mutated_defence(which_mutation_fixed);
            currentlifespan=lifespan_current;
            mutfixes(1)=mutfixes(1)+1;
        else % we move to the time point of a new bin instead
            currentbin=currentbin+1;
            newbin=1; % things will need to be recalculated in the outer loop
            currenttime=times_in_year(currentbin);
        end
    end
end
currentpopsize=popsize(currentbin);
currentbodysize=bodysizes(currentbin);
currentextmort=extrinsic(currentbin);
current_LS_extmort=LSextmort(currentbin);
LScancer=expectedhealthylifespanCancer(defence,nsteps,currentbodysize,rateonco);
currentlifespan_old=expectedhealthylifespan(defence,currentextmort,nsteps,currentbodysize,rateonco); 
currentlifespan=interp1(d,lifespan_d,defence);
current_od=optimaldefences(currentbin);
current_oFitness=optimalFitness(currentbin);
current_oLifespan=optimalLS(currentbin);
if cost==0
        costcoef=1; % if there is no cost, no multiplier for the lifespans in the next eq.
elseif cost~=0
        costcoef=(1-defence^(1/cost)); % mutiplier with the defence
end
currentfitness=costcoef*currentlifespan;
cancerp=cel(currentbodysize,nsteps,currentextmort,rateonco,defence,1000);
lagFit=currentfitness/current_oFitness;
lagD=defence-current_od;
lagD2=defence/current_od;
lagLS=currentlifespan/current_oLifespan;
DataEnd=[currenttime defence currentbodysize currentpopsize currentextmort currentlifespan currentfitness mutationrate nsteps rateonco current_LS_extmort LScancer current_od current_oFitness current_oLifespan cancerp mutfixes(1) lagFit lagLS lagD lagD2]; 
dlmwrite(filename,DataEnd,'-append'); 