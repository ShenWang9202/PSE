function x = laprnd(a,b,m,n)
% a: location parameter; b: scale parameter; m: row; n: column
% expectation: a; variance: 2b^2
% Reference: http://en.wikipedia.org/wiki/Laplace_distribution
u = rand(m,n)-0.5;
x = a - b*sign(u).*log(1-2*abs(u));
end

% cs = laprnd(40,1,1000,1)
% var(cs)
% mean(cs)