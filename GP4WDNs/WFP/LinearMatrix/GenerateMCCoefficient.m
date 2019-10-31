function [Coefficient_uncer,Variance]=GenerateMCCoefficient(OldCoefficient,MC_times,multiple)
[~,n]=size(OldCoefficient);
Coefficient_uncer=zeros(MC_times,n);
Variance = [];
for i = 1:n
    mu = OldCoefficient(i);
    sigma = multiple*mu/2.576;
    if(mu~=0)
        r = normrnd(mu,sigma,[1,MC_times]);
        Coefficient_uncer(:,i) = r';  
        Variance = [Variance; sigma^2];
    else
        Variance = [Variance; 0];
    end
end