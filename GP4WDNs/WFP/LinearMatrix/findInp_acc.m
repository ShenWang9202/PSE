function [inpname, acc] = findInp_acc(TestCase)
if(TestCase == 1)
    %inpname='Anytown2.inp';
    inpname='AnytownModify2_single.inp';
    acc = 5;
end

if(TestCase == 3)
    %inpname='tutorial5node.inp';
    inpname='ctownwithoutanyvalve.inp';
    acc = 3;
end


if(TestCase == 4)
    %inpname='tutorial5node.inp';
    inpname='ctownwithoutanyvalve_single3.inp';
    acc = 3;
end

if(TestCase == 5)
    %inpname='tutorial5node.inp';
    inpname='tutorial8node_single.inp';
    acc = 3;
end

if(TestCase == 6)
    %inpname='tutorial5node.inp';
    inpname='tutorial_lps_2PUMPS_single.inp';
    acc = 3;
end

if(TestCase == 17)
    %inpname='tutorial5node.inp';
    inpname='BAK1.inp';
    acc = 1;
end

if(TestCase == 18)
    %inpname='tutorial5node.inp';
    inpname='PES1.inp';
    acc = 1;
end

if(TestCase == 19)
    %inpname='tutorial5node.inp';
    inpname='EXN1.inp';
    acc = 3;
end

if(TestCase == 20)
    %inpname='tutorial5node.inp';
    inpname='NPCL1.inp';
    acc = 4;
end

if(TestCase == 21)
    %inpname='tutorial5node.inp';
    inpname='OBCL1.inp';
    acc = 1;
end

if(TestCase == 23)
    %inpname='tutorial5node.inp';
    inpname='Threenodes-gp_single.inp';
        acc = 2;
end

if(TestCase == 24)
    %inpname='tutorial5node.inp';
    inpname='tutorial8node5.inp';
    acc = 3;
end

if(TestCase == 25)
    %inpname='tutorial5node.inp';
    inpname='Balerma.inp';
    acc = 3;
end


