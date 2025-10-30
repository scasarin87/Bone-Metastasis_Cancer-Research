%% Define the in vivo cell line data to use as comparison
cell_line = 'renca';
follow_up = 360;

%% Compute the mean and the standard deviation

% I'll consider the average of the N simulations 
MEAN1=mean(tumor_area); %I HAVE TO CHANGE HERE THE NAMES OF THE VARIABLES I WANT TO COMPARE
STD1=std(tumor_area);


FLAG_EXPERIMENTAL_DATA = 1; % Set it to 1 if you want to compare the results with the experimental control data

%% Plot Color Matrix Definition

c = struct('rr', [0.9047, 0.1918, 0.1988], ...  %Your required color
    'bb', [0.2941, 0.5447, 0.7494], ... %Your required color
    'um', [0.0824, 0.1294, 0.4196], ... %ultra marine
    'lb', [0.1059, 0.9020, 0.9216], ... %light blue
    'lg', [0.1255, 0.9020, 0.2353], ... %light green
    'or', [1.0000, 0.5961, 0.0039], ... %ultra marine
    'ly', [0.9412, 0.9020, 0.0863], ... %ultra marine
    'br', [0.6510, 0.5725, 0.3412], ... %bronze
    'grey', [0.8863,0.8863,0.8863], ... %Your required color
    'gl', [0.8314, 0.7020, 0.7843] );   %greyed lavender

%% In vivo Experimental Data
%time=[1, 5*24, 8*24, 12*24, 15*24];

% for i = 1:size(tumor_area, 1)
%     hold on 
%     plot(1:1:follow_up, tumor_area(i, :), 'Color', c.grey, 'LineWidth', 0.08)
% end

if FLAG_EXPERIMENTAL_DATA == 1
    if strcmp(cell_line, 'c42b')
        MEDIA=[1 7.128197054 16.08757084 29.31455403 45.79737765]; % Normalized in vivo Tumor Area at 5 different time points
        STD=[0 3.719900472 5.890299748 10.73679398 16.58286703];   % STD at the same time points
        time=[1, 5*24, 8*24, 12*24, 15*24];                        % Time Points [hours]
        Interpolation = interp1(time,MEDIA,1:1:360,'pchip');       % Interpolation of the control data
        hold on
        plot(1:1:360,Interpolation,'r','LineWidth',1.2)            % Interpolated curve plot
        hold on
        errorbar(time,MEDIA,STD,'o','HandleVisibility','off')      % Errorbar definition
    elseif strcmp(cell_line, 'pc3')
        % Define PC3 Control Data and STD
        CONTROL_MEDIA = [1 8.0622 20.768 69.061 132.573]; 
        CONTROL_STD = [0 4.141 14.23 41.47 62.536];  
        % Time points
        time_points_pc3=[0, 4*24, 8*24, 17*24, 21*24];  
        % Get interpolated curves
        control_curve_pc3 = interp1(time_points_pc3, CONTROL_MEDIA, 1:1:21*24, 'pchip'); 
        % I need the curve for the first 15 days
%         control_curve_pc3 = control_curve_pc3(1:360);
%         CONTROL_MEDIA = CONTROL_MEDIA(1:4);
%         CONTROL_STD = CONTROL_STD(1:4);
%         time = time_points_pc3(1:4);
        time = time_points_pc3;
        time(1) = 1;
        hold on 
        plot(1:1:follow_up,control_curve_pc3,'b','LineWidth',1.2) 
        hold on
        errorbar(time,CONTROL_MEDIA,CONTROL_STD,'o','HandleVisibility','off')
        
    elseif strcmp(cell_line, 'renca')

        M = [1.000, 4.232, 15.115, 95.449, 309.456];
        T = [1, 96, 192, 288, 360];
        control = interp1(T, M, 1:1:360, 'pchip');

        CONTROL_MEDIA = [1, 17.612, 32.674, 141.202];
        CONTROL_STD = [0, 10.125, 22.602, 81.016];
        time_points_renca = [0, 7, 10, 15] * 24;
        control_curve_renca = interp1(time_points_renca, CONTROL_MEDIA, 1:1:15*24, 'pchip');
        time = time_points_renca;
        time(1) = 1;
        hold on 
        plot(1:1:follow_up,control_curve_renca,'r','LineWidth',1.2) 
        hold on
        errorbar(time,CONTROL_MEDIA,CONTROL_STD,'o','HandleVisibility','off')
    end
end

%% Plot of the results
% Curve 1

hold on
% 
plot(1:1:follow_up,MEAN1,'b','LineWidth',1.2)
MEAN1_grafico=zeros(1,size(time,2));
STD1_grafico=zeros(1,size(time,2));
MEAN1_grafico=MEAN1(time);
STD1_grafico=STD1(time);
hold on
errorbar(time,MEAN1_grafico,STD1_grafico,'o','HandleVisibility','off');
%ylim([0 100]);

% Curve 2
hold on
% 
% plot(1:1:360,MEAN2,'g','LineWidth',1.2)
% MEAN2_grafico=zeros(1,size(time,2));
% STD2_grafico=zeros(1,size(time,2));
% MEAN2_grafico=MEAN2(time);
% STD2_grafico=STD2(time);
% hold on
% errorbar(time,MEAN2_grafico,STD2_grafico,'o','HandleVisibility','off');
% Curve 3

% hold on
% 
% plot(1:1:360,MEAN3,'k','LineWidth',1.2)
% MEAN3_grafico=zeros(1,size(time,2));
% STD3_grafico=zeros(1,size(time,2));
% MEAN3_grafico=MEAN3(time);
% STD3_grafico=STD3(time);
% hold on
% errorbar(time,MEAN3_grafico,STD3_grafico,'o','HandleVisibility','off');

% for i = 1:size(tumor_area, 1)
%     hold on 
%     plot(1:1:follow_up, tumor_area(i, :), 'Color', c.grey, 'LineWidth', 0.1)
% end

%% Axis, Title and Legend Set up

%grid on
xlabel('Time [hours]');
ylabel('Tumor Size Normalized [ ]');
%legend('In Silico Average Control Data','Rad-223 Therapy','Cabozantinib Therapy','Combinatorial Therapy','FontSize',10)%,'Cabozantinib: size [40x24]','Cabozantinib: size [120x105]','FontSize',10)
legend('In vivo data','In silico model (avegrage trend)', '2', 'FontSize',10, 'Location', 'northwest')
% legend('In V', 'Resorption Area with ZA', 'FontSize',10)
% legend('Old Rad Algorithm', 'New Rad Algorithm', 'FontSize',10)
%legend([plot3 plot2 plot1],{'Experimental In Vivo Data','In Silico Average Control Data','In Silico Model Independent Simulations'},'FontSize',10)
title('PC3 Calibration');