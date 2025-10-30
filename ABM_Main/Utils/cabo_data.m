%% Cabozantinib Data %%

%Cabo effect on vessels taken from the Varkaris Paper 

%VASCOLARIZATION RATE WITH CABO
curvaCD31(1, :) = [104015 41283   23806   5143    4294    2468]; 
curvaCD31(2, :) = [53617	22050	20571	4566	5736	5617];
curvaCD31(3, :) = [46049	30315	24066	4547	4683	1363];
curvaCD31(4, :) = [76348	21088	21536	5249	3115	9799];

% curvaCD31(1, :) = [104015 62194   50543   38101    37535    36318]; 
% curvaCD31(2, :) = [53617	32572	31586	20916	21696	21616];
% curvaCD31(3, :) = [46049	35560	31394	18381	18476	16263];
% curvaCD31(4, :) = [76348	39508	39808	28950	27527	31983];

% curvaCD31(1, :) = [104015 72649   63910   54578    54153    53249]; 
% curvaCD31(2, :) = [53617	37834	37094	29091	29698	29635];
% curvaCD31(3, :) = [46049	38182	35058	25298	25368	23708];
% curvaCD31(4, :) = [76348	48718	48943	40799	39732	43054];

% Normalized Data wrt the initial value
c      = zeros(4, 6);
c(1, :)= curvaCD31(1, :) ./ curvaCD31(1, 1); 
c(2, :)= curvaCD31(2, :) ./ curvaCD31(2, 1);
c(3, :)= curvaCD31(3, :) ./ curvaCD31(3, 1);
c(4, :)= curvaCD31(4, :) ./ curvaCD31(4, 1);

Curvacd31 = mean(c);
timecd31  = [1,	2*24,	4*24,	6*24	9*24	12*24];
% Percentage of the vessels surviving at each hour
CURVACD31 = interp1(timecd31,Curvacd31,1:1:360,'pchip');

%Tunel is representative of apoptosis rate, si vede come aumenta con l'aumentare dell'azione del CABO
CurvaTunel=[2	2	7	110	92	176 
            6	0	7	185	52	257
            5	4	5	192	113	269
            4	2	15	85	161	332
            2	3	8	67	181	178];
        
t     = zeros(5,6);        
t(1,:)= CurvaTunel(1,:)./CurvaTunel(1,1);
t(2,:)= CurvaTunel(2,:)./CurvaTunel(2,1);
t(3,:)= CurvaTunel(3,:)./CurvaTunel(3,1);
t(4,:)= CurvaTunel(4,:)./CurvaTunel(4,1);
t(5,:)= CurvaTunel(5,:)./CurvaTunel(5,1);

Curvatunel = mean(t);
timetunel  =[1	24*2	24*4	24*6	24*9	24*21];
CURVATUNEL = interp1(timetunel,Curvatunel,1:1:24*21,'pchip'); %Polymomial Interpolation
 
time_cabo0 = [1 2*24 5*24 7*24 9*24 12*24];
 
growth_cabo_mean = [1	1.786309595	2.275766755	2.147384688	1.496097084	1.689687934]; %Normalized area in time
growth_cabo_std  = [0	0.566613971	1.229694552	0.679685014	0.749844328	0.633074577];
CURVACABO        = interp1(time_cabo0,growth_cabo_mean,1:1:24*12,'pchip');
