load('MonteCarloDemand.mat')
load('MonteCarloData.mat')
h23 = MCSolution(1,:) - MCSolution(2,:);
h37 = MCSolution(2,:) - MCSolution(6,:);
h78 = MCSolution(6,:) - MCSolution(8,:);
h_pump = MCSolution(1,:) - MCSolution(7,:);
headincreaseMC = h_pump - h23 -h37 - h78;
disp('MC')
cov(headincreaseMC,demand_MC(1,:))
cov(headincreaseMC,demand_MC(2,:))
cov(headincreaseMC,demand_MC(3,:))
cov(headincreaseMC,demand_MC(4,:))
cov(headincreaseMC,demand_MC(5,:))

k12 = A(9,9);
k23 = A(9,1);
k37 = A(9,4);
k78 = A(9,8);
b12 = Linear_pump(1,2);

headincrease_LinearTheory1 = -k12.*MCSolution(17,:) + b12 - k23.*MCSolution(9,:) - k37.*MCSolution(12,:) - k78.*MCSolution(16,:) ;
disp('LinearTheory')
cov(headincrease_LinearTheory1,demand_MC(1,:))
cov(headincrease_LinearTheory1,demand_MC(2,:))
cov(headincrease_LinearTheory1,demand_MC(3,:))
cov(headincrease_LinearTheory1,demand_MC(4,:))
cov(headincrease_LinearTheory1,demand_MC(5,:))


PumpEquation=IndexInVar.PumpEquation;
h0 =  PumpEquation(1);
r = PumpEquation(2);
nu = PumpEquation(3);
pumpsolution = 976.846008300781;
new_k12 = r * nu * (pumpsolution^(nu-1));
new_h12 = h0 + r*(pumpsolution)^nu;
new_b12 = -new_k12*pumpsolution + new_h12;




K_pipe = K_estimate(IndexInVar.PipeFlowIndex);

new_k23 = K_pipe(1);
new_k37 = K_pipe(4);
new_k78 = K_pipe(8);




headincrease_LinearTheory2 = new_k12.*MCSolution(17,:) + new_b12  - new_k23.*MCSolution(9,:) - new_k37.*MCSolution(12,:) - new_k78.*MCSolution(16,:) ;
disp('NEWLinearTheory')
cov(headincrease_LinearTheory2,demand_MC(1,:))
cov(headincrease_LinearTheory2,demand_MC(2,:))
cov(headincrease_LinearTheory2,demand_MC(3,:))
cov(headincrease_LinearTheory2,demand_MC(4,:))
cov(headincrease_LinearTheory2,demand_MC(5,:))



pumpsolution = 976.846008300781;
new1_k12 =new_k12/2;
new1_h12 = h0 + r*(pumpsolution)^nu;
new1_b12 = -new1_k12*pumpsolution + new1_h12;

headincrease_LinearTheory2 = new1_k12.*MCSolution(17,:) + new1_b12  - new_k23.*MCSolution(9,:) - new_k37.*MCSolution(12,:) - new_k78.*MCSolution(16,:) ;
disp('NEW ajust LinearTheory')
cov(headincrease_LinearTheory2,demand_MC(1,:))
cov(headincrease_LinearTheory2,demand_MC(2,:))
cov(headincrease_LinearTheory2,demand_MC(3,:))
cov(headincrease_LinearTheory2,demand_MC(4,:))
cov(headincrease_LinearTheory2,demand_MC(5,:))

h =[]
q = []
for q = 0:1:1200
    h = [h h0 + r*q^nu];
end
q = 0:1:1200;
plot(q,h)
hold on

y =new_k12.*q + new_b12;

plot(q,y)

hold on

y1 =new1_k12.*q + new1_b12;

plot(q,y1)




