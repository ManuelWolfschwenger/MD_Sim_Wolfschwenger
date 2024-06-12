clc
clear all

mu0 = 4*pi*10^(-7);
rMag = 5e-9;
Ms = 4.5e5;
N = 500;
volFrac = 0.005;

volMag = 4/3*pi*rMag^3;
magMom = volMag*Ms;
M = N*magMom^2;
Fal = 3*mu0*magMom^2/(2*pi*(2*rMag)^4);
L = (N*volMag/volFrac)^(1/3);
rc = 0.5*L;

fac = 1e-7; %mu0/(4*pi)
c = 8*M*sqrt(2*rc^3/(15*N*L^3));
d = rc^2;
delta = 1e-5;

%% real space error
Cc = @(x) 4*(x*rc).^4 + 6*(x*rc).^2 + 3;
Dc = @(x) 8*(x*rc).^6 + 20*(x*rc).^4 + 30*(x*rc).^2 + 15;
fRealExFun = @(x) fac*M./sqrt(L^3*x.^4*rc^9*N).*sqrt(13/6*Cc(x).^2 + 2/15*Dc(x).^2 - 13/15*Cc(x).*Dc(x)).*exp(-x.^2*rc^2)/Fal - delta/sqrt(2);

%% root for exact function
x = linspace(0,1e8,1e4);
Cc = 4*(x*rc).^4 + 6*(x*rc).^2 + 3;
Dc = 8*(x*rc).^6 + 20*(x*rc).^4 + 30*(x*rc).^2 + 15;
fRealEx = fac*M./sqrt(L^3*x.^4*rc^9*N).*sqrt(13/6*Cc.^2 + 2/15*Dc.^2 - 13/15*Cc.*Dc).*exp(-x.^2*rc^2)/Fal - delta/sqrt(2);

for i = 1:10000
    if fRealEx(i) > 0
        continue;
    else
        idxSmall = i;
        break;
    end
end

idxBig = idxSmall - 1;
y2 = fRealEx(idxSmall);
y1 = fRealEx(idxBig);
x2 = x(idxSmall);
x1 = x(idxBig);
k = (y2-y1)/(x2-x1);
alpha = (0 - y1)/(y2-y1)*(x2-x1) + x1;
disp(alpha*L)

figure(1)
subplot(1,2,1)
fplot(fRealExFun)
hold on
yline(0);
xline(alpha);
grid on
axis([0 1e8 -1e-4 1e-3])
xlabel('\alpha')
ylabel('\Delta F in N')
title('\alpha for given error tolerance and r_c')

%% fourier space error
c1 = 8*pi*M/L^3*sqrt(2*pi/(15*N));
d1 = (pi/(alpha*L))^2;
fRez = @(x) fac*c1*alpha*x.^1.5.*exp(-d1*x.^2)/Fal - delta/sqrt(2);
f1Rez = @(x) fac*c1*alpha*(1.5*sqrt(x).*exp(-d1*x.^2) + x.^1.5.*(-d1*2*x).*exp(-d1*x.^2))/Fal;

% find maximum value for starting value
xLin = linspace(0,30,1e3);
fLin = fac*c1*alpha*xLin.^1.5.*exp(-d1*xLin.^2)/Fal - delta/sqrt(2);

maxVal = max(fLin);
maxIdx = find(fLin == maxVal);
idx = find(fLin < 0.5*maxVal);
idx = idx(idx > maxIdx);
idx = idx(1);
kc = xLin(idx);

err = 1;
i = 1;

while err > 1e-8
    g = @(x) fRez(kc) + f1Rez(kc)*(x-kc);
    %fplot(g)
    kc = kc - fRez(kc)/f1Rez(kc);
    %xline(kc);
    err = abs(fRez(kc));
    i = i + 1;
end

subplot(1,2,2)
fplot(fRez)
hold on
yline(0);
xline(kc);
grid on
axis([0 30 -1e-4 1e-3])
xlabel('k_c')
ylabel('\Delta F in N')
title('k_c for given error tolerance and r_c')

disp(kc)