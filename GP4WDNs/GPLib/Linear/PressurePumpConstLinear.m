function [CvxExpPump,ZeorVector] = PressurePumpConstLinear(X,Pressure_Pump_Matrix_Index,X0,PumpStatus,IndexInVar)
    % variable =
    % [ head of Junction;
    %   head of reservior;
    %   head of tank;
    %   flow of pipe;
    %   flow of pump;
    %   flow of valve;
    %   speed of pump;]

    CvxExpPump = [];
    % without pump, return directly.
    if(isempty(IndexInVar.PumpFlowIndex))
        disp('No pumps in this network. Returning!');
        return
    end
    
    % getting start index of head in var
    h_start_index = min(IndexInVar.JunctionHeadIndex);
    % getting end index of head in var
    if(~isempty(IndexInVar.TankHeadIndex))
        h_end_index = max(IndexInVar.TankHeadIndex);
    else % without valves, check reservoirs
        if (~isempty(IndexInVar.ReservoirHeadIndex))
            h_end_index = max(IndexInVar.ReservoirHeadIndex);
        else % without reservoirs, only Junctions
            h_end_index = max(IndexInVar.JunctionHeadIndex);
        end
    end
    % getting start and end index of flow of pumps
    q_pump_start_index = min(IndexInVar.PumpFlowIndex);
    q_pump_end_index = max(IndexInVar.PumpFlowIndex);
    % getting start and end index of speed of pumps
    s_start_index = min(IndexInVar.PumpSpeedIndex);
    s_end_index = max(IndexInVar.PumpSpeedIndex);


    % getting GP variables
    Q = X(q_pump_start_index:q_pump_end_index);
    H = X(h_start_index:h_end_index);
    S =X(s_start_index:s_end_index);
    
    % getting initial values
    s = X0(s_start_index:s_end_index);
    q_pump = X0(q_pump_start_index:q_pump_end_index);
    
    % getting Pump equations parameter
    PumpEquation = IndexInVar.PumpEquation;
    h0_vector = PumpEquation(:,1);
    r_vector = PumpEquation(:,2);
    w_vector = PumpEquation(:,3);

    a_s = h0_vector.*s;
    [m,~] = size(Pressure_Pump_Matrix_Index);
    CvxExpPump = [];
    % PumpStatus 
    for i = 1:m
        if(PumpStatus(i))% if pump status open
            DeliveryIndex = Pressure_Pump_Matrix_Index(i,2);
            SuctionIndex = Pressure_Pump_Matrix_Index(i,1);
            %             constraints = [constraints;
            %                 H(DeliveryIndex)^(-1) * H(SuctionIndex) <= 1;% h_d >= h_s
            %                 H(DeliveryIndex) * H(SuctionIndex)^(-1) * base^(-PumpEquation(i,1)) <= 1;]; % h_d <= h_s + HeadIncreasing

            monomial_temp = 0 ;
            monomial_temp = monomial_temp + H(SuctionIndex) - H(DeliveryIndex);
            % make sure q_pump is non-negative when initializing X0 IN
            % GP_WDN_Automation.m file
            CvxExpPump = [CvxExpPump;monomial_temp + (r_vector(i)*q_pump(i)^(w_vector(i)-1)*s(i)^(2-w_vector(i))) * (Q(i)) + S(i) * (a_s(i))] ;
        end
    end
    ZeorVector=zeros(m,1);
end