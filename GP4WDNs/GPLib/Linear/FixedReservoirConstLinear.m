function [CvxExpFixed,FixedValue] = FixedReservoirConstLinear(W,X0,IndexInVar)

FixedValue = [];
CvxExpFixed = [];

SOURCE_HEAD = X0(IndexInVar.ReservoirHeadIndex);
% Tank min and max Head
TANK_HEAD = X0(IndexInVar.TankHeadIndex);

ind = 1;
for j = IndexInVar.ReservoirHeadIndex
    CvxExpFixed = [CvxExpFixed;W(j)];
    FixedValue = [FixedValue;SOURCE_HEAD(ind)];
    ind = ind + 1;
end

end