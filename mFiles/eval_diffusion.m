clc
clear all

restoredefaultpath
addpath('H:\PhD\simulation_model\MC_interact\src\build','-end');

MSDx = importdata('MSDx.txt');
MSDy = importdata('MSDy.txt');
MSDz = importdata('MSDz.txt');
velAcfx = importdata('velAcfx.txt');
velAcfy = importdata('velAcfy.txt');
velAcfz = importdata('velAcfz.txt');
deltaT = importdata('dtMean.txt');
data = importdata('data.txt');

%% theoretical diffusion coeff
mu0 = 1.25663706212e-6	;
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
Dnonint = kB*T/(6*pi*mu*rh);
%% diffusion coefficient with Einstein relation
timeOri = 10;

cutIdx = floor(length(MSDx)/timeOri)*timeOri;
MSDx = MSDx(1:cutIdx);
MSDy = MSDy(1:cutIdx);
MSDz = MSDz(1:cutIdx);

L = floor(length(MSDx)/timeOri);
MSDnx = zeros(L,1);
MSDny = zeros(L,1);
MSDnz = zeros(L,1);

for i = 0:timeOri-1
    MSDnx = MSDnx + (MSDx(i*L + 1 : (i+1)*L) - MSDx(i*L + 1));
    MSDny = MSDny + (MSDy(i*L + 1 : (i+1)*L) - MSDy(i*L + 1));
    MSDnz = MSDnz + (MSDz(i*L + 1 : (i+1)*L) - MSDz(i*L + 1));
end

MSDnx = MSDnx/timeOri;
MSDny = MSDny/timeOri;
MSDnz = MSDnz/timeOri;

time = ((0:L-1)*deltaT)';

%relate t to brown relaxation
time = time/tauB;

val1 = 0.3;
val2 = 0.99;

idx1 = floor(length(MSDnx)*val1); %stays the same for all MSD because all have the same length
idx2 = floor(length(MSDnx)*val2);

MSD = (MSDnx + MSDny + MSDnz)/3;

fitRes = createFit(time(idx1:idx2),MSDnx(idx1:idx2));
d = fitRes.p2;
k = fitRes.p1;
funx = @(x) k*x + d;
D_Einx = k/2;

fitRes = createFit(time(idx1:idx2),MSDny(idx1:idx2));
d = fitRes.p2;
k = fitRes.p1;
funy = @(x) k*x + d;
D_Einy = k/2;

fitRes = createFit(time(idx1:idx2),MSDnz(idx1:idx2));
d = fitRes.p2;
k = fitRes.p1;
funz = @(x) k*x + d;
D_Einz = k/2;

fitRes = createFit(time(idx1:idx2),MSD(idx1:idx2));
d = fitRes.p2;
k = fitRes.p1;
fun = @(x) k*x + d;
D_Ein = k/2;

figure(1)
subplot(2,2,1)
hold on
fplot(funx,[time(idx1), time(idx2)],'--','Linewidth',1.5,'Color','#0072BD')
px = plot(time,MSDnx,'Color','#0072BD');
fplot(funy,[time(idx1), time(idx2)],'--','Linewidth',1.5,'Color','#D95319')
py = plot(time,MSDny,'Color','#D95319');
fplot(funz,[time(idx1), time(idx2)],'--','Linewidth',1.5,'Color','#7E2F8E')
pz = plot(time,MSDnz,'Color','#7E2F8E');
fplot(fun,[time(idx1), time(idx2)],'--','Linewidth',1.5,'Color','#77AC30')
p = plot(time,MSD,'Color','#77AC30');
grid on
xlabel('t/\tau_B')
ylabel('MSD in m')
title('Diffusion: Einstein relation')
axis([-inf inf 0 inf])
legend([px, py, pz,p],{'MSDx','MSDy','MSDz','MSD'},'Location','northwest');

disp(D_Einx/tauB)
disp(D_Einx/tauB/Dnonint)

%% diff coeff with Green Kubo
tAcf = (0:length(velAcfx)-1)*deltaT;  %timevector for Acf

velAcf = (velAcfx + velAcfy + velAcfz)/3;
D_GKx = cumtrapz(tAcf,velAcfx);
D_GKy = cumtrapz(tAcf,velAcfy);
D_GKz = cumtrapz(tAcf,velAcfz);
D_GK = (D_GKx + D_GKy + D_GKz) / 3;

xStart = floor(0.95*length(D_GKx));

D_GK_mx = mean(D_GKx(xStart:end));
D_GK_my = mean(D_GKy(xStart:end));
D_GK_mz = mean(D_GKz(xStart:end));
D_GK_m = mean(D_GK(xStart:end));

subplot(2,2,3)
hold on
px = plot(tAcf,D_GKx,'Color','#0072BD');
py = plot(tAcf,D_GKy,'Color','#D95319');
pz = plot(tAcf,D_GKz,'Color','#7E2F8E');
p = plot(tAcf,D_GK,'Color','#77AC30','Linewidth',1.5);
%xline(tAcf(xStart),'-',{'Start measurement'});
% yline(D_Einx,'--','Color','#0072BD');
% yline(D_Einy,'--','Color','#D95319');
% yline(D_Einz,'--','Color','#7E2F8E');
%yline(D_Ein,'--','Color','#77AC30');
grid on
xlabel('time in s')
ylabel('D in m^2/s')
title('Diffusion: GK-relation')
axis([-inf inf 0 inf])
legend([px, py, pz],{'D_x','D_y','D_z'},'Location','southeast');

subplot(2,2,4)
hold on
px = plot(tAcf,velAcfx,'Color','#0072BD');
py = plot(tAcf,velAcfy,'Color','#D95319');
pz = plot(tAcf,velAcfz,'Color','#7E2F8E');
p = plot(tAcf,velAcf,'Color','#77AC30','Linewidth',1.5);
grid on
xlabel('time in s')
ylabel('velocity Acf')
title('Diffusion: GK-relation')
axis([-inf inf -inf inf])
legend([px, py, pz],{'velacf_x','velacf_y','velacf_z'},'Location','southeast');

function [fitresult, gof] = createFit(time, MSD)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time, MSD );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft);
end

