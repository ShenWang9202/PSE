function [NumNoneZero,Duncer,PDF,K]=GeneratePDF4Junction(DEMAND)
NumNoneZero = nnz(DEMAND);
sigma = 3;
percent = 0.1;
k = 4;
[m,n]=size(DEMAND);
PDF=zeros(NumNoneZero,k+1);
Duncer=zeros(NumNoneZero,k+1);
j = 0;
for i = 1:n
    i
    mu = DEMAND(i)
    if(mu~=0)
        j = j+1;
        step = 2*percent*mu/k;
        x = -mu*percent:step:mu*percent;
        x = x + mu;
        Duncer(j,:) = x;
        pdf = normpdf(x,mu,sigma)
        % size(x)
        PDF(j,:) = pdf;
    end
end
K = k+1;