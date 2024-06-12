clc
clear all

phi = 0.052;
Ms = 4.5e5;
temp = 300;
kB = 1.38e-23;
mu0 = 4*pi*10^(-7);

rm = @(lambda) (lambda*18*kB*temp/(mu0*Ms^2*pi))^(1/3);

%% analytical solution
vol = @(dm) pi/6*dm.^3;
lambda = @(dm) mu0*Ms^2*vol(dm)/(24*kB*temp);
lambda_iv = @(dm) (vol(dm)*Ms)^2/(dm^3*kB*temp);
lambda_wang = @(dm) (vol(dm)*Ms)^2/(4*pi*mu0*kB*temp*dm^3);

z_field = @(dm) exp(2*lambda(dm))./(3*lambda(dm).^2);
z_nofield = @(dm) exp(2*lambda(dm))./(3*lambda(dm).^3);

n_field = @(dm) 0.5 + sqrt(0.25 + phi*z_field(dm));
n_nofield = @(dm) 0.5 + sqrt(0.25 + phi*z_nofield(dm));

%rosensweig
% n_field_ros = @(dm) (1 - 2*phi*exp(2*lambda(dm))./(3*lambda(dm).^2)).^(-1);
%n_nofield_ros = @(dm) (1 - 2*phi*exp(2*lambda(dm))./(3*lambda(dm).^3)).^(-1);

%% simulation data

coords = readmatrix('coords5.txt');
data = readmatrix('data5.txt');

rMag = data(1,2);
lBox = data(1,3);
N = data(2,3); % Number of particles

vol = pi/6*(2*rMag)^3;
m = vol*Ms;

threshold = 1.5*kB*temp*lambda(2*rMag(1)); %*lambda(2*rMag(1))
disp(lambda(2*rMag))

% Initialize variables
steps = length(coords)/N;

tt = 1;

for t = 0:100:steps-1
    
    chain = zeros(N,1); % Vector to store the chain
    belong = zeros(N,1); % Vector to mark to which chain particles belong
    idx = t * N;
    c = 1;
    
    time(tt) = coords(idx + 1,13);
    
    % Loop over all particles
    for i = 1:N-1
        for j = (i + 1) : N
            
            flag = true;
            
            %interaction vector
            rij = coords(j+idx,1:3) - coords(i+idx,1:3);
            
            %periodic boundary condition
            idx_r = rij > lBox/2;
            rij = rij - idx_r*lBox;
            idx_r = rij < -lBox/2;
            rij = rij + idx_r*lBox;
            
            r = norm(rij);
            rij = rij/r;
            
            En_dd = mu0*m^2/(4*pi*r.^3)*(3*dot(coords(i+idx,7:9),rij)*dot(coords(j+idx,7:9),rij) - dot(coords(i+idx,7:9),coords(j+idx,7:9)));
            
            if En_dd >= threshold
                if  belong(i) == 0 && belong(j) == 0 && flag
                    chain(c) = 2;  % create new chain
                    belong(i) = c; % now particles belong to chain nr. c
                    belong(j) = c;
                    c = c + 1;
                    flag = false;
                end
                
                if xor(belong(i) ~= 0, belong(j) ~= 0) && flag
                    
                    %check what particle already belongs to a chain
                    if belong(i) ~= 0
                        c_idx = belong(i);
                    else
                        c_idx = belong(j);
                    end
                    
                    chain(c_idx) = chain(c_idx) + 1;
                    belong(i) = c_idx;
                    belong(j) = c_idx;
                    
                end
            end
        end
    end
    
    % add singe particles to chain
    for i = 1:N
        if belong(i) == 0
            belong(i) = c;
            chain(c) = 1;
            c = c + 1;
        end
    end
    
    chain = chain(chain > 0);
    n_chain(tt) = mean(chain);
    disp(t/steps*100)
    
    if sum(chain) ~= N 
        disp('error incorrect number of particles')
    end
    
    tt = tt + 1;
end

%% plots
figure(1)
plot(time,n_chain)
hold on
yline(n_nofield(2*rMag));
%xline(time(xStart),'--r',{'beginning of measurement', 'range'});
axis([0,inf,1,inf])
grid on
xlabel('time in s')
ylabel('mean chain length')
title('chain size over time')

% xMax = 16e-9;
% lambda_min = 2;
% xMin = (lambda_min*24*kB*temp*6/(mu0*Ms^2*pi))^(1/3);
% figure(2)
% fplot(n_field,[xMin,xMax]);
% hold on
% fplot(n_nofield,[xMin,xMax]);
% grid on
% axis([xMin,xMax,1,inf]);
% title('chain size distribution')
% xlabel('magnetic diameter in m')
% ylabel('mean chain size')
% legend('B_{mag} = 1T','B_{mag} = 0T','Location','northwest')

