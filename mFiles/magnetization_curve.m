clc
clear all

restoredefaultpath
addpath('H:\PhD\simulation_model\MC_interact\src\build','-end');

t = importdata('t.txt');
mz = importdata('mz.txt');
nz = importdata('nz.txt');

figure(1)
plot(t,mz)
hold on
plot(t,nz)
xlabel('time in s')
ylabel('magnetization')
grid on
axis([0 inf -0.1 1])