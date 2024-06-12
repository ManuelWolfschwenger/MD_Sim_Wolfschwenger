clc
clear all

restoredefaultpath
addpath('H:\PhD\simulation_model\MC_interact\src\build','-end');

data = importdata('data.txt');
Pxy = importdata('Pxy.txt');
Pxz = importdata('Pxz.txt');
Pyz = importdata('Pyz.txt');
deltaT = importdata('dtMean.txt');

idxCut = 1; %wait till equilibrium

lBox = data(1,3);          % box length
numbPart = data(2,3);      % number of atoms
T = data(5,3);             %temperature

kB = 1.38064852e-23;    %bolzmann constant
Vbox = lBox^3;

%% brown relaxation time
mu0 = 1.25663706212e-6;
rh = data(1,1);
rm = data(1,2);
T = data(5,3);
kB = 1.380648e-23;
D_anal = kB*T/(6*pi*853.8e-6*rh);
Vm = 4/3*pi*rm^3;
Ms = data(4,3);
mu = 2.414e-5*10^(247.8/(T-140));
Vh = 4/3*pi*rh^3;

lambda = mu0*(Ms*Vm)^2/(4*pi*kB*T*8*rm^3)*(rm/rh)^3;
tauB = 3*mu*Vh/(kB*T);
%% smoothing and cutting
Pxy = Pxy(idxCut:end);
Pxz = Pxz(idxCut:end);
Pyz = Pyz(idxCut:end);

Pxy = smooth(Pxy);
Pxz = smooth(Pxz);
Pyz = smooth(Pyz);

%% viscosity with Einstein relation
timeOri = 200; %number of time origins to compute average of

L = floor(length(Pxy)/timeOri); %length of series
tVisEin = (0:L-1)'*deltaT;  %timevector for Acf

%relate t to brown relaxation
tVisEin = tVisEin/tauB;

MSDxy = zeros(L,1);
MSDxz = zeros(L,1);
MSDyz = zeros(L,1);

for i = 1:timeOri
    MSDxy = MSDxy + cumtrapz(tVisEin,Pxy((i-1)*L + 1:i*L)).^2;
    MSDxz = MSDxz + cumtrapz(tVisEin,Pxz((i-1)*L + 1:i*L)).^2;
    MSDyz = MSDyz + cumtrapz(tVisEin,Pyz((i-1)*L + 1:i*L)).^2;
end

MSDxy = MSDxy/timeOri;
MSDxz = MSDxz/timeOri;
MSDyz = MSDyz/timeOri;

idx1 = floor(0.8*L);

fitRes = createFit(tVisEin(idx1:end),MSDxy(idx1:end));
d = fitRes.p2;
k = fitRes.p1;
funxy = @(x) k*x + d;
visEinxy = Vbox/(2*kB*T)*k;

fitRes = createFit(tVisEin(idx1:end),MSDxz(idx1:end));
d = fitRes.p2;
k = fitRes.p1;
funxz = @(x) k*x + d;
visEinxz = Vbox/(2*kB*T)*k;

fitRes = createFit(tVisEin(idx1:end),MSDyz(idx1:end));
d = fitRes.p2;
k = fitRes.p1;
funyz = @(x) k*x + d;
visEinyz = Vbox/(2*kB*T)*k;

% %stuff for slope plot
% step = 20;
% 
% t = tVisEin(1:step:end);
% MSDxys = MSDxy(1:step:end);
% MSDxzs = MSDxz(1:step:end);
% MSDyzs = MSDyz(1:step:end);
% kxy = Vbox/(2*kB*T)*diff(MSDxys)./diff(t);
% kxz = Vbox/(2*kB*T)*diff(MSDxzs)./diff(t);
% kyz = Vbox/(2*kB*T)*diff(MSDyzs)./diff(t);

figure(1)
subplot(2,2,1)
hold on
pxy = plot(tVisEin,MSDxy,'Color','#0072BD');
fplot(funxy,[tVisEin(idx1), tVisEin(end)],'--','Linewidth',1.5,'Color','#0072BD')
pxz = plot(tVisEin,MSDxz,'Color','#D95319');
fplot(funxz,[tVisEin(idx1), tVisEin(end)],'--','Linewidth',1.5,'Color','#D95319')
pyz = plot(tVisEin,MSDyz,'Color','#7E2F8E');
fplot(funyz,[tVisEin(idx1), tVisEin(end)],'--','Linewidth',1.5,'Color','#7E2F8E')
grid on
xlabel('t/\tau_B')
ylabel('')
legend([pxy, pxz, pyz],{'MSDxy','MSDxz','MSDyz'},'Location','northwest');
title('Viscosity: Einstein relation')
axis([0 inf -inf inf])

% subplot(2,2,2)
% hold on
% pxy = plot(t(1:end-1),kxy,'Color','#0072BD');
% yline(visEinxy,'--','Color','#0072BD');
% pxz = plot(t(1:end-1),kxz,'Color','#D95319');
% yline(visEinxz,'--','Color','#D95319');
% pyz = plot(t(1:end-1),kyz,'Color','#7E2F8E');
% yline(visEinyz,'--','Color','#7E2F8E');
% grid on
% xlabel('t/\tau_B')
% ylabel('viscosity in Pa \cdot s')
% title('Viscosity: Einstein relation')
% axis([0 inf -inf inf])
% legend([pxy, pxz, pyz],{'\eta_{xy}','\eta_{xz}','\eta_{yz}'},'Location','northwest');

%% GK viscosity
timeOri = 200; %number of timeorigins to compute average of

L = floor(length(Pxy)/timeOri); %length of series
tVisGK = (0:L-1)'*deltaT;  %timevector for Acf

pAcfxy = autocorr(Pxy,L-1)*var(Pxy);
pAcfxz = autocorr(Pxz,L-1)*var(Pxz);
pAcfyz = autocorr(Pyz,L-1)*var(Pyz);

intAcfxy = Vbox/(kB*T)*cumtrapz(tVisGK,pAcfxy);
intAcfxz = Vbox/(kB*T)*cumtrapz(tVisGK,pAcfxz);
intAcfyz = Vbox/(kB*T)*cumtrapz(tVisGK,pAcfyz);

idx1 = floor(0.8*length(tVisGK));

%relate t to brown relaxation
tVisGK = tVisGK/tauB;

subplot(2,2,3)
hold on
pxy = plot(tVisGK,intAcfxy,'Color','#0072BD');
pxz = plot(tVisGK,intAcfxz,'Color','#D95319');
pyz = plot(tVisGK,intAcfyz,'Color','#7E2F8E');
yline(visEinxy*tauB,'--','Color','#0072BD');
yline(visEinxz*tauB,'--','Color','#D95319');
yline(visEinyz*tauB,'--','Color','#7E2F8E');
xline(tVisGK(idx1),'-',{'Start measurement'});
grid on
xlabel('t/\tau_B')
ylabel('vis in Pa*s')
legend([pxy, pxz, pyz],{'\eta_{xy}','\eta_{xz}','\eta_{yz}'},'Location','northwest');
title('Viscosity: GK relation')
axis([0 inf -inf inf])

subplot(2,2,4)
hold on
pxy = plot(tVisGK,pAcfxy,'Color','#0072BD');
pxz = plot(tVisGK,pAcfxz,'Color','#D95319');
pyz = plot(tVisGK,pAcfyz,'Color','#7E2F8E');
grid on
xlabel('t/\tau_B')
ylabel('pressure tensor Acf')
title('Pressure Tensor Acf- GK relation')
axis([0 inf -inf inf])
legend([pxy, pxz, pyz],{'pacf_{xy}','pacf_{xz}','pacf_{yz}'},'Location','northeast');
%set(gca,'XScale','log');

function [fitresult, gof] = createFit(time, MSD)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time, MSD );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft);
end

