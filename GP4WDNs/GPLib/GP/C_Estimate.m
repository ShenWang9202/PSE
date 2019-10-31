function c_estimate = C_Estimate(d,X0,base,q_pipe_start_index,q_pipe_end_index)

% get data from epanet
%     h_start_index = 1;
%     h_end_index = 8;
%     q_pipe_start_index = 9;
%     q_pipe_end_index = 16;
%     q_pump_start_index = 17;
%     q_pump_end_index = 17;
%     s_start_index = 18;
%     s_end_index = 18;
%X0




q_pipe = X0(q_pipe_start_index:q_pipe_end_index,1);
FlowUnits = d.getFlowUnits;
if(strcmp('LPS',FlowUnits{1})) % convert to gpm
    L_pipe = d.getLinkLength *Constants4WDN.m2feet; % ft
    D_pipe = d.getLinkDiameter *Constants4WDN.mm2inch; % inches ; be careful, pump's diameter is 0
end
if(strcmp('GPM',FlowUnits{1})) % convert to gpm
    L_pipe = d.getLinkLength ; % ft
    D_pipe = d.getLinkDiameter; % inches ; be careful, pump's diameter is 0
end
C_pipe = d.getLinkRoughnessCoeff; % roughness of pipe
% if(strcmp('LPS',FlowUnits{1}))
%     diameter_conversion = mm2m;
%     Volum_conversion = L2cube_m;
% end

diameter_conversion = Constants4WDN.feet2inch;
Volum_conversion = Constants4WDN.GPM2CFS;


PipeIndex = d.getLinkPipeIndex;
L_pipe = L_pipe(PipeIndex);
D_pipe = D_pipe(PipeIndex)/diameter_conversion;
C_pipe = C_pipe(PipeIndex);

Headloss_pipe_R = 4.727 * L_pipe./((C_pipe*Volum_conversion).^(1.852))./(D_pipe.^(4.871));

% get the R cofficient for pipe 1 to 8
% if(strcmp('LPS',FlowUnits{1}))
%     Headloss_pipe_R = 10.66 * L_pipe./((C_pipe*Volum_conversion).^(1.852))./(D_pipe.^(4.871));
% end
c_estimate = [];
[m,n] = size(q_pipe)
for i = 1:m
    c_estimate = [c_estimate;base^((Headloss_pipe_R(i)*abs(q_pipe(i))^(0.852)-1)*q_pipe(i))];
    %c_estimate = [c_estimate;base^((Headloss_pipe_R(i)*abs(q_pipe(i))^(0.852)))];
end
end
