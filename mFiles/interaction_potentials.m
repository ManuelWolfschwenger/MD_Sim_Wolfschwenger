
clear all

%% physical constants
kB = 1.38e-23;
mu0 = 4*pi*10^(-7);
NA = 6.02214076e23;
el = 1.6e-19;
epsilon0 = 8.8541e-12;

%% assembly properties
dm = 2*4e-09;
rm = dm/2;
temp = 300;
V = 4/3*pi*rm^3;
Ms = 4.5e5;
m = Ms*V;

%% steric repulsion
zeta = 5e17;
delta = 2e-9;
d_delta = dm + 2*delta;
rd = d_delta/2;

pot_ster = @(r) (0 <= r & r <= d_delta).*pi*dm^2*zeta*kB*temp/(2*delta).*(d_delta - r.*(log(d_delta./r) + 1))/(kB*temp) + (r > d_delta).*0;
f_ster = @(r) (dm <= r & r <= d_delta).*kB*temp*pi*dm^2*zeta/(2*delta).*log(d_delta./r) + (r > d_delta).*0;
fun1 = @(r) 1 - exp(-pot_ster(r));
dh = integral(fun1,0,d_delta);

%% vdw
hamaker = 33e-21;

s = @(r) r - 2*rm;
pot_vdw = @(r) -hamaker/6*(2*rm*rm./(s(r).*(r+2*rm)) + 2*rm*rm./(r.^2 - (rm - rm)^2) + log((r.^2 - (rm+rm)^2)./(r.^2 - (rm-rm)^2)))/(kB*temp);
f_vdw = @(r) -hamaker.*r*(4*rm*rm)^3./(6*(s(r).*(r+rm+rm)).^2.*(r.^2-(rm-rm)^2).^2);

%% dip-dip
pot_dip = @(r) -2*mu0*m*m./(4*pi*r.^3)/(kB*temp);
f_dip = @(r) 3*mu0./(4*pi*r.^4)*(-2*m^2);

%% electrostat
epsilonr = 78.5;
z = 1; %valency of ions
gamma0 = 0.04; %surface potential in V
c = 0.05*1000; %concentration from mol/l to mol/m^3

epsilon = epsilon0*epsilonr;
gamma = tanh(z*el*gamma0/(4*kB*temp));
kappa = sqrt(el^2*NA*c*z^2/(epsilon*kB*temp));

pot_el = @(r) 32*pi*epsilon*rm*(kB*temp/(z*el))^2*gamma^2*exp(-kappa*s(r))/(kB*temp);
f_el = @(r) kappa*32*pi*epsilon*rd*(kB*temp/(z*el))^2*gamma^2*exp(-kappa*s(r));

%% sum
pot_sum = @(r) pot_ster(r) + pot_vdw(r) + pot_dip(r) + pot_el(r);
f_sum = @(r) f_ster(r) + f_vdw(r) + f_dip(r) + f_el(r);

%% cut off calculation
delta = 0.01;
c = 128*hamaker*rm^6*pi/(18*mu0*m^2);
fCutOff = @(r) c*r/((r^2-4*rm^2)^2) - delta; 
fCutOff1 = @(r) c*(4*rm^2+3*r^2)/(4*rm^2-r^2)^3; %derivation of cutoff

rCut = 1.1*dm;
err = 1;

%% plots
rEndForce = 1.25*d_delta;
rEndPot = 2*d_delta;

figure(2)
subplot(1,2,1)
title('interaction potential')
fplot(pot_ster,[dm,d_delta])
hold on
fplot(pot_vdw,[dm,rEndPot])
fplot(pot_dip,[dm,rEndPot])
fplot(pot_el,[dm,rEndPot])
fplot(pot_sum,[dm,rEndPot],'Linewidth',1.5)
xl1 = xline(d_delta,'--',{'shell surface'},'Linewidth',1);
xl2 = xline(dm,'--',{'core surface'},'Linewidth',1);
xl1.LabelVerticalAlignment = 'top';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'top';
xl2.LabelHorizontalAlignment = 'left';

grid on
legend('steric','vdW','dip-dip','el.-stat','sum','Location','northeast')
title('potential energy')
xlabel('distance between particle centers in m')
ylabel('interaction energy in k_BT units')
axis([0.95*dm rEndPot -30 50])

subplot(1,2,2)
fplot(f_ster,[dm,rEndForce])
hold on
fplot(f_vdw,[dm,rEndForce])
fplot(f_dip,[dm,rEndForce])
fplot(f_el,[dm,rEndForce])
fplot(f_sum,[dm,rEndForce],'Linewidth',2)
xl1 = xline(d_delta,'--',{'shell surface'},'Linewidth',1);
xl2 = xline(dm,'--',{'core surface'},'Linewidth',1);
xl1.LabelVerticalAlignment = 'top';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment = 'top';
xl2.LabelHorizontalAlignment = 'left';
grid on
legend('ster','vdw','dip','el','sum','Location','northeast')
title('interaction forces')
xlabel('distance between particle centers in m')
ylabel('interaction force in N')
axis([0.95*dm rEndForce -1e-10 2.5e-10])

% figure(2) 
% fplot(fCutOff,[dm,5*dm], 'LineWidth', 3);
% hold on
% 
% while err > 1e-8
%     g = @(r) fCutOff(rCut) + fCutOff1(rCut)*(r-rCut);
%     fplot(g)
%     rCut = rCut - fCutOff(rCut)/fCutOff1(rCut);
%     xline(rCut);
%     err = abs(fCutOff(rCut));
% end
% 
% xl1 = xline(d_delta,'--',{'shell surface'},'Linewidth',1);
% xl2 = xline(dm,'--',{'core surface'},'Linewidth',1);
% xl3 = xline(rCut,'--',{'cut off'},'Linewidth',1);
% xl1.LabelVerticalAlignment = 'top';
% xl1.LabelHorizontalAlignment = 'left';
% xl2.LabelVerticalAlignment = 'top';
% xl2.LabelHorizontalAlignment = 'left';
% xl3.LabelVerticalAlignment = 'top';
% xl3.LabelHorizontalAlignment = 'left';
% grid on
% axis([0.95*dm 5*dm -0.1 0.2])
% xlabel('distance of centerpoints in m')
% ylabel('F_{vdw}/F_{dip}')
% title('cut off for short range potential')
% 
% disp(rCut)


