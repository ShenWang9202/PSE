function [CvxExpHeadLoss,CEstimateVetor] = PressurePipeConstLinear(X,Pressure_Pipe_Matrix_Index,X0,d,IndexInVar)
    % variable =
    % [ head of Junction;
    %   head of reservior;
    %   head of tank;
    %   flow of pipe;
    %   flow of pump;
    %   flow of valve;
    %   speed of pump;]

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
    
    % getting start and end index of flow of pipes
    q_pipe_start_index = min(IndexInVar.PipeFlowIndex);
    q_pipe_end_index = max(IndexInVar.PipeFlowIndex);
    % getting GP variables
    Q = X(q_pipe_start_index:q_pipe_end_index);
    H = X(h_start_index:h_end_index);
    % getting  C_Estimate
    CEstimateVetor = CEstimateLinear(d,X0,q_pipe_start_index,q_pipe_end_index);
    [m,~] = size(Pressure_Pipe_Matrix_Index);

    CvxExpHeadLoss = [];
    for i = 1:m
        monomial_temp = 0;
        monomial_temp = monomial_temp + H(Pressure_Pipe_Matrix_Index(i,1)) - H(Pressure_Pipe_Matrix_Index(i,2));
        %CvxExpHeadLoss = [CvxExpHeadLoss;monomial_temp - (Q(i)) - c_estimate_value(i)] ;
        CvxExpHeadLoss = [CvxExpHeadLoss;monomial_temp - Q(i)] ;
    end
end