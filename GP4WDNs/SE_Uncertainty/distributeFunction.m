close all;
mu = 200;
sigma = 5;
figure
x = mu + (-20:1:20)
y1 = normpdf(x,mu,sigma)
plot(x,y1)
hold on

x = mu+50+ (-20:1:20)
y2 = normpdf(x,mu+50,sigma)
plot(x,y2)

y3 = conv(y1,y2);
figure
plot(y3);

pd = makedist('Normal','mu',mu,'sigma',sigma);

y_cdf = cdf(pd,x);
hold on
plot(x,y_cdf)

likelihood = unidpdf(5,4:9)
plot(likelihood)


x = normrnd(61.7,4.487927836,1000,1);
y = normrnd(75,7.5,1000,1);
z = normrnd(82.5,5.04,1000,1);
sumR = x-y-z;
mean(sumR), sqrt(var(sumR))


Demo_dMATLAC.m

