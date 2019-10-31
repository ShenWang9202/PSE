function [AnalysisMatrix,B] = Construct_Variance_A_b(NumberofX,A,demand_MC)

xx = sym('x',[1 NumberofX]);


randValue  = randperm(NumberofX);
sqrtValue = sqrt([1:double(NumberofX)]);
xValue = randValue.*sqrtValue;

expMatrix = [];
for i = 1:NumberofX
    Ai = A(i,:);
    expMatrix=[expMatrix xx*Ai'];
end

[m,n_d]=size(demand_MC);
PesudoDemand = zeros(NumberofX,n_d);
for i = 1:m
    PesudoDemand(i,:) = demand_MC(i,:);
end
B = [];
expression = [];
combination =[];

for i = 1:NumberofX
    for j = i:NumberofX
        cov_value = cov(PesudoDemand(i,:),PesudoDemand(j,:));
        B = [B;cov_value(1,2)];
    end
end
tic
for i = 1:NumberofX
    for j = i:NumberofX
        combination = [combination xValue(i)*xValue(j)] ;
    end
end
toc
tic
for i = 1:NumberofX
    for j = i:NumberofX
        expression = [expression expMatrix(i)*expMatrix(j)] ;
    end
end
toc
% for i = 1:NumberofX
%     for j = i:NumberofX
%         combination = [combination xx(i)*xx(j)] ;
%         expression = [expression expMatrix(i)*expMatrix(j)] ;
%         cov_value = cov(PesudoDemand(i,:),PesudoDemand(j,:));
%         B = [B;cov_value(1,2)];
%     end
% end

n = (NumberofX+1)*NumberofX*0.5;
AnalysisMatrix = zeros(n,n);

val=243; %value index of which to be found
[uar,ia,ib]=unique(combination);
sind=1; lind=length(uar);
while lind-sind>0,
mind=ceil((sind+lind)/2);
if val>=uar(mind), sind=mind; else, lind=mind-1; end
end
ind=ia(lind);


[~,n] = size(expression);
for i = 1:n
    tic
    result = expression(i);
    [cxy, txy] = coeffs(result,xx);
    txyValue = subs(txy,xx,xValue)
    [~,m] = size(cxy);
    toc
    tic
    for j = 1:m
        ind = find(combination == txyValue(j));
 %       AnalysisMatrix(i,ind) = cxy(j);
    end
    toc
    i
end
% 
% A= [-1 1 0 0 0 0 0 0 1];
% e1 =  xx * A'
% result = e1*e1;
% expand(result)
% [cxy, txy] = coeffs(result,xx)
% 
% find(combination == txy(4))

end