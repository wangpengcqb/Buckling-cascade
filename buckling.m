
% xx = Untitled(:,6);
% yy = Untitled(:,7);
% zz = Untitled(:,8) - Untitled(:,5);
% cftool(xx,yy,zz)
% 
% x = -100:0.1:100;
% y = 0:0.1:90;
% 
% [X,Y] = meshgrid(x,y);
% Z = fittedmodel(X,Y);

for i = 1:length(Z)
   
    ztemp = Z(:,i);
    clear u3 x2;
    nc = 1;
    for j=1:length(ztemp)  
        if isfinite(ztemp(j))
            u3(nc,1) = ztemp(j);
            x2(nc,1) = y(j);
            nc = nc + 1;
        end
    end
if(nc>1)    
nor = length(u3);
u3_f = u3;

width = x2(end) - x2(1);

dw(1) = (u3(2) - u3(1))/(x2(2) - x2(1));
for m = 2:nor-1
    dw(m) = (u3(m+1) - u3(m-1))/(x2(m+1) - x2(m-1));    
end
dw(nor) = (u3(nor) - u3(nor-1))/(x2(nor) - x2(nor-1)); 

widthc = 0;
for m = 1:nor-1
    widthc = widthc + 1/2*(sqrt(1+dw(m)^2) + sqrt(1+dw(m+1)^2))*(x2(m+1)-x2(m));
end

eps(i,1) = 1 - width/widthc;

for m = 1:nor-1
u3_f(nor+m) = - u3_f(nor-m);
end
u3_f(2*nor-1) = [];

n = 2*nor-2;
fs = 1/(y(2)-y(1));

f = (0:n/2-1)/n*fs;
f = f';

yf = fft(u3_f);
TT = abs(yf(1:n/2));

ave_f = sqrt(norm(TT.*f)^2/norm(TT)^2);
ave_lm(i,1) = 1./ave_f;

A_ave(i,1) = norm(TT)*2/n;

% resu = [TT' f'];
% 
% amp = [TT(2) TT(4) TT(6) TT(8) TT(10) TT(12) TT(14) TT(16) TT(18)];

% hold off;
% plot(f,TT)
% saveas(gcf,'fourier')  
else
    eps(i,1) = 0;
    ave_lm(i,1) = 0;
    A_ave(i,1) = 0;;
end
       
end