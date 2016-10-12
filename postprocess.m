clear all;

load x2.txt;
load u3.txt;

x2(:,1) = [];
u3(:,1) = [];


for i = 1:21
    hold on;
 plot(x2(i,:),u3(i,:))
 
end

saveas(gcf,'profile')

nor = length(u3(21,:));

u3_f = u3(21,:);

for i = 1:nor-1
u3_f(nor+i) = - u3_f(nor-i);
end
u3_f(2*nor-1) = [];

n = 2*nor-2;
fs = 1/0.01;

f = (0:n/2-1)/n*fs;

yy = fft(u3_f);
TT = abs(yy(1:n/2));

ave_f = sqrt(norm(TT.*f)^2/norm(TT)^2);
ave_lm = 1./ave_f;

resu = [TT' f'];

amp = [TT(2) TT(4) TT(6) TT(8) TT(10) TT(12)];

hold off;
plot(f,TT)
saveas(gcf,'fourier')