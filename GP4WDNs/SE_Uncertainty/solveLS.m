% 
% A=[1 2 -2 3; 2 4 -3 4; 5 10 -8 11];
% b=[2 5 12]';
% 
% A =[
%     1 -1 0 0 0 0 0 0 0;
%     0 1 -1 0 -1 0 0 0 0;
%     0 0 1 -1 0 1 0 0 0;
%     0 0 0 1 0 0 0 -1 0;
%     0 0 0 0 0 1 1 1 0;
%     0 0 0 0 1 0 -1 0 -1;];
% b = [100; 100; 200; 20; 390; 28];
% 
% format rat;
% [S_H, S_P]=solveLS(A,b)



function [S_H, S_P] = solveLS(A,b)
if size(A,1) ~= length(b) 
    error('Parameter error');
    return;
else
    B = [A,b];  
    rank_A = rank(A);   
    rank_B = rank(B);   
    if rank_A ~= rank_B 
        disp('no solution');
        S_H = [];
        S_P = [];
    else if rank_B == size(A,2)  
            %size(A,2) 
            disp('unique solution');
            S_P = A\b;  
            S_H = [];
        else
            disp('infinite solution');
            S_H = null(A,'r'); 
            S_P = A\b; 
        end
    end
end
