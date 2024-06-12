clc
clear all

restoredefaultpath
addpath('H:\PhD\simulation_model\QuantMRX\src\build','-end');

t = importdata('t.txt');
mz = importdata('mz.txt');

figure(1)
plot(t(1:50:end),mz(1:50:end))
xlabel('time in s')
ylabel('magnetization')
grid on