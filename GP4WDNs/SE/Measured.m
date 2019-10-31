function [A_over,b_over] = Measured(A,b,MeasuredIndex,MeasuredValue)
    MeasuredIndex = MeasuredIndex(:);
    MeasuredValue = MeasuredValue(:);
    [m,n] = size(MeasuredIndex);
    [m2,n2] = size(MeasuredValue);
    if(m~=m2 && n~=n2)
        disp('something is wrong!');
        return;
    end

    [~,nA] = size(A);
    A_added = zeros(n,nA);
    b_added = zeros(n,1);
    for i = 1:n
        A_added(i,MeasuredIndex(i)) = 1;
    end
    b_added = MeasuredValue;
    A_over = [A;A_added];
    b_over = [b;b_added];
end