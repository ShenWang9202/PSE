function [constraints,c_estimate,ValveFlagOut] = PressureValveConst(X,Pressure_Valves_Matrix_Index,X0,base,NodeElevations,ValveSetting,ValveTypesString,verify,ValveStatus,ValveFlag,Solution,FlowUnits,iter,IndexInVar)
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

    q_valve_start_index = min(IndexInVar.ValveFlowIndex);
    q_valve_end_index = max(IndexInVar.ValveFlowIndex);
    %     q_pipe_start_index = min(IndexInVar.PipeFlowIndex);
    %     q_pipe_end_index = max(IndexInVar.PipeFlowIndex);
    %     s_start_index = min(IndexInVar.PumpSpeedIndex);
    %     s_end_index = max(IndexInVar.PumpSpeedIndex);
    
    Q = X(q_valve_start_index:q_valve_end_index);
    H = X(h_start_index:h_end_index);
    q_valve = X0(q_valve_start_index:q_valve_end_index);
    head_s = Solution(h_start_index:h_end_index,:);
    if(strcmp('GPM',FlowUnits))
        conversion = Constants4WDN.PSIperFT;
        M2FT = 1;
    end

    if(strcmp('LPS',FlowUnits))
        conversion = 1;
        M2FT = Constants4WDN.M2FT;
    end
    
    [~,n] = size(ValveTypesString);
    c_estimate = [];
    constraints = [];
    ValveFlagOut = ValveFlag;
    for i = 1:n
        if (strcmp('PRV',ValveTypesString{i}))
            if(ValveStatus(i)==0) %  the PRV is off.
                constraints = [constraints;
                    Q(i) == base^(0);  % set q_{valve} = 0;
                    ];
            end
            if(ValveStatus(i))%  the PRV is working (on or active).
                DownStreamIndex = Pressure_Valves_Matrix_Index(i,2);
                UpStreamIndex = Pressure_Valves_Matrix_Index(i,1);
                elevation = NodeElevations(DownStreamIndex);
                constraints = [constraints;
                    Q(i)^(-1) <= 1; % flow q_{valve} >= 0
                    ];
                if(ValveFlag(i) == 2)
                    constraints = [constraints;H(DownStreamIndex) == base^((elevation + ValveSetting(i)/conversion)*M2FT)];
                end
                if(ValveFlag(i) == 1)
                    constraints = [constraints;H(DownStreamIndex) == H(UpStreamIndex)];
                    % The first two conditions are for stablizing ; the
                    % last one is when the pressure is greater than PRV
                    % settings, then set it.
                    
                    % this is for MPC not for Water flow problem
%                     if( iter>=2 && abs(head_s(DownStreamIndex,end)-head_s(DownStreamIndex,end-1)) < 1 && head_s(DownStreamIndex,end)/M2FT - elevation > ValveSetting(i)/conversion)
%                         ValveFlagOut(i) = 2;
%                     end
                end
                
            end
        end

        if (strcmp('FCV',ValveTypesString{i}))
            if(ValveStatus(i)) % open
                DownStreamIndex = Pressure_Valves_Matrix_Index(i,2);
                UpStreamIndex = Pressure_Valves_Matrix_Index(i,1);
                constraints = [constraints;H(DownStreamIndex) == H(UpStreamIndex)];
                %constraints = [constraints;Q(i) * base^(-FCVSetting) <= 1]; % flow always >=0
            else
                constraints = [constraints;Q(i) == base^(0);]; % flow ==0
            end

        end

        if (strcmp('GPV',ValveTypesString{i}))
            c_estimate_value = C_Estimate_Valve(base,q_valve(i));
            monomial_temp = H(Pressure_Valves_Matrix_Index(i,1))^(1) * H(Pressure_Valves_Matrix_Index(i,2))^(-1);
            f_mono = monomial_temp * (Q(i))^(-1) * c_estimate_value^(-1);
            constraints = [constraints;f_mono == 1];
            c_estimate = [c_estimate;c_estimate_value];
        end
        %         if (strcmp('PSV',ValveTypesString{i}))
        %             ;
        %         end
        %         if (strcmp('PBV',ValveTypesString{i}))
        %             ;
        %         end
        %         if (strcmp('TCV',ValveTypesString{i}))
        %             ;
        %         end
    end

end