clc
clear all

coords = readmatrix('coords.txt');
data = readmatrix('data.txt');

rHydr = data(:,1);
rMag = data(:,2);
lBox = data(1,3);
numbPart = data(2,3);

time = coords(1:numbPart:end,13);
maxPairs = numbPart * (numbPart - 1) / 2;

% evaluation of steps over time
steps = length(coords)/numbPart; %steps stored from c++

% properties of histogram
bins = 500;
rangeRdf = 8*rMag(1);
tMin = 2000;%2000;
tMax = 3000;%3000;

cols = tMax - tMin + 1;
binWidth = rangeRdf/bins;
histData = zeros(bins,cols);
normFac = lBox^3/(2*pi*binWidth^3*numbPart^2);

o = 1;
j = 1;
for t = tMin : tMax
    for j1 = 1 : numbPart
        for j2 = j1 + 1 : numbPart - 1
            rij = coords(j1 + t*numbPart,1:3) - coords(j2 + t*numbPart,1:3);
            
            %periodic boundary condition
            idx_r = rij > lBox/2;
            rij = rij - idx_r*lBox;
            idx_r = rij < -lBox/2;
            rij = rij + idx_r*lBox;
            
            rij = norm(rij);
            
            if rij < 2*rMag(1)
                overlap(o) = abs(rij-2*rMag(1))/(2*rMag(1));
                o = o + 1;
                rij = 2*rMag(1);
            end
            
            if rij < rangeRdf
                n = floor(rij / binWidth) + 1;
                histData(n,j) = histData(n,j) + 1;
            end
        end
    end
    
    for k = 1:bins
        histData(k,j) =  histData(k,j)*normFac/((k - 0.5)^2);
    end
    
    disp(j)
    j = j+1; 
end

histData = mean(histData,2);
x = ((linspace(1,bins,bins) - 0.5)*binWidth)/(2*rMag(1))';
%disp(max(overlap))

figure(1)
plot(x,histData)
hold on
grid on
yline(1,'--');
xlabel('distance in d_m')
ylabel('radial distribution function')
axis([0, inf, 0, 10])

%% plot hist data of all three
data = xlsread('histdata.xlsx','A3:D502');
x = data(:,1);
y_ster = data(:,2);
y_el = data(:,3);
y_elster = data(:,4);

figure(1)
hold on
plot(x,y_ster,'Color','#0072BD')
plot(x,y_el,'Color','#D95319')
plot(x,y_elster,'Color','#7E2F8E')
grid on
yline(1,'--');
xlabel('distance in d_m')
ylabel('radial distribution function g(r)')
axis([0, 4, 0, 10])
legend('steric','electrostat.','electrosteric')


