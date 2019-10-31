function constraints = FlowRangeConst(X,X0,verify,z,IndexInVar)
    % variable =
    % [ head of Junction;
    %   head of reservior;
    %   head of tank;
    %   flow of pipe;
    %   flow of pump;
    %   flow of valve;
    %   speed of pump;]

    % For nodes like source or tanks shouldn't have mass equations.
    % getting start index of flow in var
    q_start_index = min(IndexInVar.PipeFlowIndex);
    % getting end index of flow in var
    if(~isempty(IndexInVar.ValveFlowIndex))
        q_end_index = max(IndexInVar.ValveFlowIndex);
    else % without valves
        if (~isempty(IndexInVar.PumpFlowIndex))
            q_end_index = max(IndexInVar.PumpFlowIndex);
        else % without pumps
            q_end_index = max(IndexInVar.PipeFlowIndex);
        end
    end
    % getting GP variables
    Q = X(q_start_index:q_end_index);
    q_pipe = X0(q_start_index:q_end_index,1);
    [m,~] = size(q_pipe);
    f_mono = [];
    monomial_temp = 1.0;
    for i = 1:m
        monomial_temp = monomial_temp * Q(i)^(-q_pipe(i));
    end
    if(verify == 0)
        constraints = monomial_temp <= z^(-1);
    else
        constraints = f_mono;
    end

end