A_test = [1 -2 1;
    1 -1 0;
    0 1 -1;
    1 0 0;
    0 0 1];
b_test = [100;0;0;0;0];
w = [0.01 0.01 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1;]
A_test'*w*A_test\(A_test'*w*b_test)

A_test\b_test

