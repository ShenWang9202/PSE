clc
clear
TestCase = 23;
[inpname, acc] = findInp_acc(TestCase);
[d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution] = PrepareA(inpname,TestCase);


Headloss_pipe_R = PipeCoeff(ForConstructb);
R23 = Headloss_pipe_R(1);
if(~isempty(IndexInVar.PumpEquation))
    h0_vector = IndexInVar.PumpEquation(:,1);
    r_vector = IndexInVar.PumpEquation(:,2);
    w_vector = IndexInVar.PumpEquation(:,3);
end
% initial value
q12 = 0;
q23 = 0;

% water demand at Junction 2
d2 = 100;

% measurement at all nodes
h1 = 700;
h2 = 910;
h3 = 908;
%
delta_h31= h3-h1;
delta_h21= h2-h1;
delta_h23= h2-h3;

%% sufficient scenario

%objective = @(x) (R23*x(1)*abs(x(1))^(0.852) + R34*x(2)*abs(x(2))^(0.852) - delta_h)^2;
objective = @(x) (-R23*x(2)*abs(x(2))^(0.852) + h0_vector + r_vector*(x(1))^(w_vector) - delta_h31)^2;

% initial guess
x0 = [q12;q23];

% variable bounds
lb = [0 -1000];%-1000.0 * ones(2);
ub = [1000 1000];%1000.0 * ones(2);

% show initial objective
disp(['Initial Objective: ' num2str(objective(x0))])

% linear constraints
A = [];
b = [];
Aeq = [1 -1];
beq = [d2];

% nonlinear constraints
nonlincon = [];

% optimize with fmincon
%[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] 
% = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon);

% show final objective
disp(['Final Objective: ' num2str(objective(x))])

% print solution
disp('Sufficient Scenerio Solution')
disp(['q12 = ' num2str(x(1))])
disp(['q23 = ' num2str(x(2))])


%% over-determinded scenario
 
%objective = @(x) ((-R23*x(2)*abs(x(2))^(0.852) + delta_h23)^2 +(h0_vector + r_vector*(x(1))^(w_vector) - delta_h21)^2);
objective = @(x) ( (x(1)-h2)^2 + (x(2)-h1)^2 + (x(3)-h3)^2 + (-R23*x(5)*abs(x(5))^(0.852) + h2-h3)^2 +(h0_vector + r_vector*(x(4))^(w_vector) - (x(1)-x(2)))^2);
% initial guess
x0 = [h2;h1;h3;q12;q23];

% variable bounds
lb = [0 0 0 0 -1000];%-1000.0 * ones(2);
ub = [1000 1000 1000 1000 1000];%1000.0 * ones(2);

% show initial objective
disp(['Initial Objective: ' num2str(objective(x0))])

% linear constraints
A = [];
b = [];
Aeq = [0 0 0 1 -1];
beq = [d2];

% nonlinear constraints
nonlincon = [];

% optimize with fmincon
%[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] 
% = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon);

% show final objective
disp(['Final Objective: ' num2str(objective(x))])
% print solution
disp('0ver-determined Scenerio Solution')
disp(['h2 = ' num2str(x(1))])
disp(['q12 = ' num2str(x(4))])
disp(['q23 = ' num2str(x(5))])
h0_vector + r_vector*(x(1))^(w_vector)