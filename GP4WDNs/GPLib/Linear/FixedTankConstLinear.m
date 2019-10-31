function [CvxExpFixed,FixedValue] = FixedTankConstLinear(W,X0,IndexInVar)

FixedValue = [];
CvxExpFixed = [];
% Tank min and max Head
TANK_HEAD = X0(IndexInVar.TankHeadIndex);

% fix the head of tank
ind = 1;
for j = IndexInVar.TankHeadIndex
    CvxExpFixed = [CvxExpFixed;W(j)];
    FixedValue = [FixedValue;TANK_HEAD(ind)];
    ind = ind +1;
end


end