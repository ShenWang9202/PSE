function  [d,InitialParameter,SettingsNStatus,IndexInVar,MassEnergyMatrix4GP,Variable_Symbol_Table,Solution]=Prepare(inpname,TestCase)
d=epanet(inpname);
% d.plot('nodes','yes','links','yes','highlightnode',{'1','8'},'highlightlink',{'7'},'fontsize',8);

Velocity=[];
Pressure=[];
T=[];
Demand=[];
Head=[];
Flows=[];
TankVolume=[];
HeadLoss = [];
LinkSettings = [];
LinkStatus = [];
% Another way to Simulate all
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1;
d.getTimeHydraulicStep
while (tstep>0)
    t=d.runHydraulicAnalysis;   %current simulation clock time in seconds.
    Velocity=[Velocity; d.getLinkVelocity];
    Pressure=[Pressure; d.getNodePressure];
    Demand=[Demand; d.getNodeActualDemand];
    TankVolume=[TankVolume; d.getNodeTankVolume];
    HeadLoss=[HeadLoss; d.getLinkHeadloss];
    Head=[Head; d.getNodeHydaulicHead];
    Flows=[Flows; d.getLinkFlows];
    T=[T; t];
    LinkSettings = [LinkSettings;d.getLinkSettings];
    LinkStatus = [LinkStatus;d.getLinkStatus];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis

%% Get Solution from EPANET
[m,n] = size(T);

%PipeIndex = d.getLinkPipeIndex; This method would miss the CVPIPE, e.g
%P446 in CTOWN
% So there is a bug in getLinkPipeIndex and getLinkPipeNameID functions,
% these two funcitons missed CVPIPE type.

PipeIndex = 1:d.getLinkPipeCount;
PumpIndex = d.getLinkPumpIndex;
ValveIndex = d.getLinkValveIndex;
NodeJunctionIndex = d.getNodeJunctionIndex;
ReservoirIndex = d.getNodeReservoirIndex;
NodeTankIndex = d.getNodeTankIndex;

% Pump and Valve Status
PumpStatus = LinkStatus(:,PumpIndex);
ValveStatus = LinkStatus(:,ValveIndex);

% Settings for all types of links
PipeRoughness = LinkSettings(:,PipeIndex);
PumpSpeed = LinkSettings(:,PumpIndex); % take effect when pumpstatus is open
ValveSettings = LinkSettings(:,ValveIndex);

PumpSpeedNew = PumpSpeed .* PumpStatus;
%result =xor(PumpSpeedNew,PumpStatus)


% Generate Solution for later validation.
Solution = [];
for i = 1:m
    Solution = [Solution; [Head(i,:) Flows(i,:) PumpSpeedNew(i,:)]];
end
Solution = Solution';
%% Generate Demand for NodeJuncion;
Demand_known = [];
saved_temp_index = [];
Demand_NodeJunction = Demand(:,NodeJunctionIndex);
[m,n] = size(Demand_NodeJunction);
Demand_known = [Demand_known;Demand_NodeJunction(1,:)];
saved_temp_index = [saved_temp_index 1];
distance_demand = 0;
for i = 2:m
    distance_demand = norm(Demand_NodeJunction(i,:) - Demand_known(end,:));
    if(distance_demand >= 0)
        Demand_known = [Demand_known;Demand_NodeJunction(i,:)];
        saved_temp_index = [saved_temp_index i];
    end
end
Demand_known = Demand_known';
% find corresponding Solution
Solution = Solution(:,saved_temp_index);
ValveStatus = ValveStatus';
ValveStatus = ValveStatus(:,saved_temp_index);
PumpStatus = PumpStatus';
PumpStatus = PumpStatus(:,saved_temp_index);


%% Generate Mass and Energy Matrice
NodeNameID = d.getNodeNameID; % the Name of each node   head of each node
LinkNameID = d.getLinkNameID; % the Name of each pipe   flow of each pipe

NodesConnectingLinksID = d.getNodesConnectingLinksID; %
[m,n] = size(NodesConnectingLinksID);
NodesConnectingLinksIndex = zeros(m,n);

for i = 1:m
    for j = 1:n
        NodesConnectingLinksIndex(i,j) = find(strcmp(NodeNameID,NodesConnectingLinksID{i,j}));
    end
end
NodesConnectingLinksIndex
% Generate MassEnergyMatrix
[m1,n1] = size(NodeNameID);
[m2,n2] = size(LinkNameID);
MassEnergyMatrix = zeros(n2,n1);

for i = 1:m
    MassEnergyMatrix(i,NodesConnectingLinksIndex(i,1)) = -1;
    MassEnergyMatrix(i,NodesConnectingLinksIndex(i,2))= 1;
end
% Display
MassEnergyMatrix


%% Generate Mass Matrix
% For nodes like source or tanks shouldn't have mass equations.
MassMatrix = MassEnergyMatrix(:,NodeJunctionIndex)';
[m,~] = size(MassMatrix)
MassMatrixIndexCell = cell(m,2);
[RowPos,ColPos] = find(MassMatrix == 1);
[m,~] = size(RowPos);
for i = 1:m
    MassMatrixIndexCell(RowPos(i),1) = {[MassMatrixIndexCell{RowPos(i),1},ColPos(i)]};
end

[RowNeg,ColNeg] = find(MassMatrix == -1);
[m,~] = size(RowNeg);
for i = 1:m
    MassMatrixIndexCell(RowNeg(i),2) = {[MassMatrixIndexCell{RowNeg(i),2},ColNeg(i)]};
end



%% Generate Energy Matrix

% for pipe
EnergyPipeMatrix = -MassEnergyMatrix(PipeIndex,:);

[m,~] = size(EnergyPipeMatrix);
EnergyPipeMatrixIndex = zeros(m,2);
[RowPos,ColPos] = find(EnergyPipeMatrix == 1);
[m,~] = size(RowPos);
for i = 1:m
    EnergyPipeMatrixIndex(RowPos(i),1) = ColPos(i);
end

[RowNeg,ColNeg] = find(EnergyPipeMatrix == -1);
[m,~] = size(RowNeg);
for i = 1:m
    EnergyPipeMatrixIndex(RowNeg(i),2) = ColNeg(i);
end

% for Pump
EnergyPumpMatrix = -MassEnergyMatrix(PumpIndex,:);

[m,~] = size(EnergyPumpMatrix);
EnergyPumpMatrixIndex = zeros(m,2);
[RowPos,ColPos] = find(EnergyPumpMatrix == 1);
[m,~] = size(RowPos);
for i = 1:m
    EnergyPumpMatrixIndex(RowPos(i),1) = ColPos(i);
end

[RowNeg,ColNeg] = find(EnergyPumpMatrix == -1);
[m,~] = size(RowNeg);
for i = 1:m
    EnergyPumpMatrixIndex(RowNeg(i),2) = ColNeg(i);
end

% for valve
EnergyValveMatrix = -MassEnergyMatrix(ValveIndex,:);


[m,~] = size(EnergyValveMatrix);
EnergyValveMatrixIndex = zeros(m,2);
[RowPos,ColPos] = find(EnergyValveMatrix == 1);
[m,~] = size(RowPos);
for i = 1:m
    EnergyValveMatrixIndex(RowPos(i),1) = ColPos(i);
end

[RowNeg,ColNeg] = find(EnergyValveMatrix == -1);
[m,n] = size(RowNeg);
for i = 1:m
    EnergyValveMatrixIndex(RowNeg(i),2) = ColNeg(i);
end

MassEnergyMatrix4GP = struct('MassMatrixIndexCell',{MassMatrixIndexCell},...
    'EnergyPipeMatrixIndex',EnergyPipeMatrixIndex,'EnergyPumpMatrixIndex',EnergyPumpMatrixIndex,...
    'EnergyValveMatrixIndex',EnergyValveMatrixIndex);

%% Find the index for variable.
% variable =
% [ head of Junction;
%   head of reservior;
%   head of tank;
%   flow of pipe;
%   flow of pump;
%   flow of valve;
%   speed of pump;]

% count of each element
PipeCount = d.getLinkPipeCount;
PumpCount = d.getLinkPumpCount;
ValveCount = d.getLinkValveCount;

JunctionCount = d.getNodeJunctionCount;
ReservoirCount = d.getNodeReservoirCount;
TankCount = d.getNodeTankCount;

% total count.
NumberofX = PipeCount + PumpCount + ValveCount;
NumberofX = NumberofX + JunctionCount + ReservoirCount + TankCount;
NumberofX = NumberofX + PumpCount; % speed of pump

% index for each element.
JunctionHeadIndex = 1:JunctionCount;

BaseCount4Next = JunctionCount;
ReservoirHeadIndex = (BaseCount4Next+1):(BaseCount4Next+ReservoirCount);

BaseCount4Next = BaseCount4Next+ReservoirCount;
TankHeadIndex = (BaseCount4Next+1):(BaseCount4Next + TankCount);

BaseCount4Next = BaseCount4Next+TankCount;
PipeFlowIndex = (BaseCount4Next+1):(BaseCount4Next + PipeCount);

BaseCount4Next = BaseCount4Next+PipeCount;
PumpFlowIndex = (BaseCount4Next+1):(BaseCount4Next + PumpCount);

BaseCount4Next = BaseCount4Next+PumpCount;
ValveFlowIndex = (BaseCount4Next+1):(BaseCount4Next + ValveCount);

BaseCount4Next = BaseCount4Next+ValveCount;
PumpSpeedIndex = (BaseCount4Next+1):(BaseCount4Next + PumpCount);


%% Calculate accurate pumpequation coefficiencies

% The paramter of pump curve from EPANET, BUT the  r_vector is not accurate, we
% need to calibrate them via the solution provide by EPANET.
if(TestCase == 1) %AnytownModify
    PumpEquation = [300 -3.9549e-06 1.91;]; %3.9549e-06 -4.059E-06
end
if(TestCase == 2) %BWSN_TestCase_1
    PumpEquation = [445.00 -1.947E-05 2.28;
        740.00 -8.382E-05 1.94;
        ];
end
if(TestCase == 3) %ctown
    PumpEquation = [229.659 -0.005969 1.36;
        229.659 -0.005969 1.36;
        295.28 -0.0001146 2.15;
        393.7 -3.746E-006 2.59;
        295.28 -0.0001146 2.15;
        295.28 -4.652E-005 2.41;
        ];
end
if(TestCase == 4 )
    PumpEquation = [393.7008 -3.8253E-006 2.59;];
    %PumpEquation = [200*0.3048 -0.01064 2;];
end

if(TestCase == 5 )
    PumpEquation = [393.7008 -3.8253E-006 2.59;];
    %PumpEquation = [200*0.3048 -0.01064 2;];
end
if(TestCase == 6 )
    PumpEquation = [200 -5.952E-6 2;];
    %PumpEquation = [200*0.3048 -0.01064 2;];
end
if(TestCase ==7  )
    PumpEquation = [393.7008 -3.8253E-006 2.59;];
    %PumpEquation = [200*0.3048 -0.01064 2;];
end

if(TestCase ==19  )
    PumpEquation = [393.7008 -3.8253E-006 2.59;];
    %PumpEquation = [200*0.3048 -0.01064 2;];
end
% PumpEquation
% find the head index
%     variable =
%     [ head of Junction;
%       head of reservior;
%       head of tank;
%       flow of pipe;
%       flow of pump;
%       flow of valve;
%       speed of pump;]
% getting start index of head in var
h_start_index = min(JunctionHeadIndex);
% getting end index of head in var
if(~isempty(TankHeadIndex))
    h_end_index = max(TankHeadIndex);
else % without valves, check reservoirs
    if (~isempty(ReservoirHeadIndex))
        h_end_index = max(ReservoirHeadIndex);
    else % without reservoirs, only Junctions
        h_end_index = max(JunctionHeadIndex);
    end
end
% getting start and end index of flow of pumps
q_pump_start_index = min(PumpFlowIndex);
q_pump_end_index = max(PumpFlowIndex);
% getting start and end index of speed of pumps
s_start_index = min(PumpSpeedIndex);
s_end_index = max(PumpSpeedIndex);

FlowUnits = d.getFlowUnits;
if(strcmp('LPS',FlowUnits{1})) % convert to gpm
   head_unit_conversion = Constants4WDN.M2FT;
   flow_unit_conversion = Constants4WDN.LPS2GMP;
end
if(strcmp('GPM',FlowUnits{1})) % convert to gpm
    head_unit_conversion = 1;
    flow_unit_conversion = 1;
end

HeadofAllNode = Solution(h_start_index:h_end_index,1)*head_unit_conversion;
SpeedofAllPump = Solution(s_start_index:s_end_index,1);
FlowofAllPump = Solution(q_pump_start_index:q_pump_end_index,1)*flow_unit_conversion;

h0_vector = PumpEquation(:,1);
w_vector = PumpEquation(:,3);

[m,~] = size(EnergyPumpMatrixIndex);
% search for each pump  and update  r_vector
for i = 1:m
%     get the index of junction node connecting pumps.
    DeliveryIndex = EnergyPumpMatrixIndex(i,2);
    SuctionIndex = EnergyPumpMatrixIndex(i,1);
    HeadIncreaseofPump= HeadofAllNode(DeliveryIndex) - HeadofAllNode(SuctionIndex);
    
    r_vector = (h0_vector(i)-HeadIncreaseofPump/(SpeedofAllPump(i)*SpeedofAllPump(i))) * (SpeedofAllPump(i)/FlowofAllPump(i))^(w_vector(i));
    PumpEquation(i,2)= -r_vector;
end

%%
IndexInVar = struct('NumberofX',NumberofX,'JunctionCount',JunctionCount,'JunctionHeadIndex',JunctionHeadIndex,...
    'ReservoirHeadIndex',ReservoirHeadIndex,...
    'TankHeadIndex',TankHeadIndex,'PipeFlowIndex',PipeFlowIndex,'PumpFlowIndex',PumpFlowIndex,...
    'ValveFlowIndex',ValveFlowIndex,'PumpSpeedIndex',PumpSpeedIndex,...
    'PumpEquation',PumpEquation);

%% Variable_Symbol_Table


Variable_Symbol_Table = cell(NumberofX,2);
temp_i = 1;
NodeIndexInVar = d.getNodeIndex;
LinkIndexInVar = d.getLinkIndex  + d.getNodeCount;
LinkPumpNameID = d.getLinkPumpNameID;
for i = NodeIndexInVar
    Variable_Symbol_Table{i,1} =  NodeNameID{temp_i};
    temp_i = temp_i + 1;
end

temp_i = 1;
for i = LinkIndexInVar
    Variable_Symbol_Table{i,1} =  LinkNameID{temp_i};
    temp_i = temp_i + 1;
end
% remove speed, not a optimization variable any more.
temp_i = 1;
for i = PumpSpeedIndex
    Variable_Symbol_Table{i,1} = strcat('Speed_',LinkPumpNameID{temp_i});
    temp_i = temp_i + 1;
end

for i = 1:NumberofX
    Variable_Symbol_Table{i,2} =  strcat('W_',int2str(i));
end


%% Initial parameter 
%Elevation
Elevation = d.getNodeElevations;
PipeFLowAverage = mean(Solution,2);
ReservoirHead = Head(:,ReservoirHeadIndex);
TankHead = Head(:,TankHeadIndex);
%%
%FlowUnits = {'GPM'};
FlowUnits = d.getFlowUnits;
if(strcmp('LPS',FlowUnits{1}))
    M2FT = Constants4WDN.M2FT;
    LPS2GMP = Constants4WDN.LPS2GMP;
end

if(strcmp('GPM',FlowUnits{1}))
    M2FT = 1;
    LPS2GMP = 1;
end
Demand_known = Demand_known * LPS2GMP;




InitialParameter = struct('Elevation',Elevation,'PipeFLowAverage',PipeFLowAverage,...
    'ReservoirHead',ReservoirHead,'Demand_known',Demand_known,'M2FT',M2FT,'LPS2GMP',LPS2GMP);

SettingsNStatus = struct('PipeRoughness',PipeRoughness,...
    'PumpStatus',PumpStatus,'PumpSpeed',PumpSpeed,...
    'ValveSettings',ValveSettings,...
    'ValveStatus',ValveStatus,...
    'FlowUnits',FlowUnits);
% pump status(open or close) should be viewed as known, and shouldn't be
% placed in Solution, but it is possible that the status of Pumps can be
% variables, so we just viewed them as variables with fixed value now.

