clear
s = rng;
RandNum = [];
for i = 1:1000
r = normrnd(3,1,[1,5])
RandNum = [RandNum;r];
end
%RandNum=sort(RandNum);
figure
histogram(RandNum(:,1),20);
hold on
histogram(RandNum(:,2),20);
hold on
histogram(RandNum(:,3),20);


Fs = 3000; % sampling rate of 1000 Hz
x = RandNum(:,1);
xdft = fft(x);
xdft = xdft(1:floor(length(x)/2)+1);
DF = Fs/length(x); % frequency increment
freqvec = 0:DF:Fs/2;
figure
plot(freqvec,abs(xdft))