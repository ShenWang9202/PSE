function [D_uncer,DemandIndex,Variance] = GenerateMCDemand(demand,MC_times,HeadIndex,multiple,flowConverter)
[~,n]=size(demand);
D_uncer=zeros(MC_times,n);
DemandIndex = [];
Variance = [];
for i = 1:n
    mu = demand(i);
    sigma = multiple*mu/2.576;
    if(mu~=0)
        r = normrnd(mu,sigma,[1,MC_times]);
        D_uncer(:,i) = r';  
        DemandIndex = [DemandIndex HeadIndex(i)];
        sigma = sigma *flowConverter; % flowConverter is for LPS inp files. demand will convert later.
        Variance = [Variance; sigma^2];
    else
        Variance = [Variance; 0];
    end
end
