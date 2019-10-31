function [Measurement_uncer,Variance]=GenerateMCError(measurments,MC_times,accuracy)
[~,n]=size(measurments);
Measurement_uncer=zeros(MC_times,n);
Variance = [];
% the accuracy of sensor are 2%
if (accuracy == 0)
    return;
end
for i = 1:n
    mu = measurments(i);
    sigma = accuracy*mu/2.576;
    if(mu~=0)
        r = normrnd(mu,sigma,[1,MC_times]);
        Measurement_uncer(:,i) = r';  
        Variance = [Variance; sigma^2];
    else
        Variance = [Variance; 0];
    end
end
%Measurement_uncer = Measurement_uncer';