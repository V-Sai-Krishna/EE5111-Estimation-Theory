pd = makedist('tLocationScale','mu',0.5,'sigma',1,'nu',1);
t = truncate(pd,0,1);
x=linspace(0,2,1000);
% plot(x,pdf(t,x));
pdist=random(t,1,N);
hpdist = histcounts(pdist, hgrid, 'Normalization', 'probability');
            cdfhpdist = cumsum(hpdist);
            plot(cdfhpdist,[1:1:1001]);