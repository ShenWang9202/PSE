function c_estimate = C_Estimate_Valve(base,q_valve)
%FlowUnits = d.getFlowUnits;
Headloss_valve_R = 4.127E-4;
c_estimate = base^((Headloss_valve_R*abs(q_valve)^(0.852)-1)*q_valve);
end
