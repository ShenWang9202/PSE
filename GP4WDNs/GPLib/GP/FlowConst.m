function constraints = FlowConst(X,MassMatrixCell,Demand,verify,IndexInVar)
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
    [m,~] = size(MassMatrixCell);
    f_mono = [];
    for i = 1:m
        monomial_temp = 1.0;
        PostiveMatrixIndex = cell2mat(MassMatrixCell(i,1));
        [~,n] = size(PostiveMatrixIndex);
        for j = 1:n
            monomial_temp = monomial_temp * (Q(PostiveMatrixIndex(j)))^(1);
        end
        NegativeMatrixIndex = cell2mat(MassMatrixCell(i,2));
        [~,n] = size(NegativeMatrixIndex);
        for j = 1:n
            monomial_temp = monomial_temp * (Q(NegativeMatrixIndex(j)))^(-1);
        end
        f_mono = [f_mono; Demand(i)^(-1)* monomial_temp] ;
    end
    [m,~] = size(f_mono);
    if(verify == 0)
        constraints = f_mono == ones(m,1);
    else
        constraints = f_mono;
    end

end