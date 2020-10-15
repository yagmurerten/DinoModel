% code for simulating different 'lag' scenarios for Fig 4 in the main text

logcosts=linspace(-4.1,-0.1,100);
costs=10.^logcosts;

k=1e-4;n=3;
% we look at different reductions in extrinsic mortality
extred=linspace(1,50,100);
% using either the ancestral dinosaur size or the final bird size
ancestralbs=220.7;newbs=0.8;

ancestralextmort=getExtMort(ancestralbs);newextmort=getExtMort(newbs);
scaled_new_e=newextmort*ones(1,100)./extred;

for i=1:100
    c=costs(i);
    [ancestrald,fitness1,lifespan1]=getoptimalLSdefence(ancestralextmort,n,k,ancestralbs,c);
    [newd,fitness2,lifespan2]=getoptimalLSdefence(newextmort,n,k,newbs,c);
    % calculating the lags between
    % (1) bird lineage ending with the ancestral defences
    % (2) bird lineage ending with pre-flight shrunk defences
    % and the optimal level of defences for the given extrinsic mortality
    % and body size
    for j=1:100
        e=scaled_new_e(j);
        [postredd,postredfit,postredLS]=getoptimalLSdefence(e,n,k,newbs,c);
        L1=expectedhealthylifespan(ancestrald,e,n,newbs,k);
        L2=expectedhealthylifespan(newd,e,n,newbs,k);
        costcoef1=1-ancestrald.^(1/c); 
        costcoef2=1-newd.^(1/c); 
        fit1=costcoef1.*L1; fit2=costcoef2.*L2;
        lagd1=ancestrald-postredd;
        lagd2=newd-postredd;
        lagfit1=fit1/postredfit;
        lagfit2=fit2/postredfit;
        lagLS1=L1/postredLS;
        lagLS2=L2/postredLS;
        temp=horzcat(c,extred(j),lagd1,lagd2,lagfit1,lagfit2,lagLS1,lagLS2);
        dlmwrite('lags.txt',temp,'-append');
    end
end

