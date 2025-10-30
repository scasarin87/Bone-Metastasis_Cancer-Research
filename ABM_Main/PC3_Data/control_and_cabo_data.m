% Define PC3 Control Data and STD
CONTROL_MEDIA = [1 6.579 23.689 49.612 61.934 113.800]; 
CONTROL_STD=[0 2.041 13.725 26.117 34.407 43.536];  

% Define PC3 Cabo Data and STD
CABO_MEDIA = [1 5.381 17.809 30.836 59.570 123.146]; 
CABO_STD=[0 1.676 5.241 15.839 42.903 71.834];

% Time points
time=[0, 4*24, 8*24, 11*24, 17*24, 21*24];  

% Get interpolated curves
control_curve = interp1(time, CONTROL_MEDIA, 1:1:21*24, 'pchip');    
cabo_curve = interp1(time, CABO_MEDIA, 1:1:21*24, 'pchip');    

% Plot curve
plot(1:1:15*24,control_curve(1:360),'r','LineWidth',1.2)            
hold on
plot(1:1:15*24,cabo_curve(1:360),'b','LineWidth',1.2)            
hold on
errorbar(time(1:4),CONTROL_MEDIA(1:4),CONTROL_STD(1:4),'o','HandleVisibility','off')     
hold on
errorbar(time(1:4),CABO_MEDIA(1:4),CABO_STD(1:4),'o','HandleVisibility','off') 

xlabel('Time [h]');
ylabel('Tumor Size Normalized [ ]');
legend('PC3 Control','PC3 Cabo', 'FontSize',10)
