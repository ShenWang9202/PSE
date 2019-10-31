clc;
clear;
close all;
TestCase = 1;

if(TestCase == 1)
    inpname='tutorial8node.inp';
end

if(TestCase == 5)
    %inpname='tutorial5node.inp';
    inpname='tutorial8node.inp';
end

if(TestCase == 7)
   inpname='Threenodes-gp.inp';
   %inpname='tutorial4price19_copy.inp';
end

 [d,IndexInVar,InitialParameter,ForConstructA,ForConstructb,Variable_Symbol_Table,Solution,MassEnergyMatrix4GP]=PrepareA(inpname,TestCase)

 
 PumpEquation =IndexInVar.PumpEquation;
plotpump(PumpEquation)


Headloss_pipe_R = PipeCoeff(ForConstructb);
plotpipe(Headloss_pipe_R)
