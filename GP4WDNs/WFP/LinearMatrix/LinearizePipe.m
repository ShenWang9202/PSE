function [K_pipe K_c] = LinearizePipe(forConstructb,Pipeflow)
FlowUnits = forConstructb.FlowUnits;
if(strcmp('LPS',FlowUnits)) % convert to gpm
    L_pipe = forConstructb.LinkLength *Constants4WDN.m2feet; % ft
    D_pipe = forConstructb.LinkDiameter *Constants4WDN.mm2inch; % inches ; be careful, pump's diameter is 0
end
if(strcmp('GPM',FlowUnits)) % convert to gpm
    L_pipe = forConstructb.LinkLength ; % ft
    D_pipe = forConstructb.LinkDiameter; % inches ; be careful, pump's diameter is 0
end
C_pipe = forConstructb.LinkRoughnessCoeff; % roughness of pipe

diameter_conversion = Constants4WDN.feet2inch;
Volum_conversion = Constants4WDN.GPM2CFS;

PipeIndex = forConstructb.LinkPipeIndex;
L_pipe = L_pipe(PipeIndex);
D_pipe = D_pipe(PipeIndex)/diameter_conversion;
C_pipe = C_pipe(PipeIndex);


%Headloss_pipe_R = 4.727 * L_pipe./((C_pipe*Volum_conversion).^(1.852))./(D_pipe.^(4.871));
common = 4.727 * L_pipe./((Volum_conversion).^(1.852))./(D_pipe.^(4.871));
q = Pipeflow';
K_pipe =1.852.*common.*C_pipe.^(-1.852).*(abs(q).^(0.852));
K_c = -1.852.*common.*(q.^1.852).*C_pipe.^(-2.852);
% get the R cofficient for pipe 1 to 8
% if(strcmp('LPS',FlowUnits{1}))
%     Headloss_pipe_R = 10.66 * L_pipe./((C_pipe*Volum_conversion).^(1.852))./(D_pipe.^(4.871));
% end
end
