%% Six Hours Angiogenesis %%

%  Here we perform the most recurrent vessel generation, each 6 hours new
%  vessels are generated according to the current tumor area

% Compute tumor area at the current iteration
curr_tumor_area = sum(sum(bone == site.tumor_edge | bone == site.tumor | bone == site.vessel | bone == site.vessel_cabo)) * site_dim * site_dim; %Area in squared µm           
% Compute the vessels to be generated with this angiogenetic step
new_vess_number = round((Mean_Vessel_Density) * curr_tumor_area) - vessels_number; 

% Vessels mask creation: finds the core of pre-established vessels;
% generates new ones outside of it
mask = ones(rows, columns); 
vessels_mask_creation; % individuates a core of already established vessels. i wil try to generate new vessels outside it

% I wanto to compute the tumor baricenter
stats = regionprops(bone == site.tumor_edge | bone == site.tumor | bone == site.vessel | bone == site.vessel_cabo); %Statistics of the tumor geometry
centroids = cat(1, stats.Centroid); %Computation of the tumor centroids coordinates
r_c = fix(centroids(1, 2));
c_c = fix(centroids(1, 1)); 

% Creation of eigth part tumor masks
maschere_ott = octants_mask_creation(r_c, c_c, rows, columns); % Creation of eigth part tumor masks
mask_center_vessel  = zeros(rows, columns); %It will contain the points being  vessels centers

% I compute for each eighth part the cells increment wrt 5
% hours before. The probability of creating a new vessel in
% that region is linearly dependent on this increment wrt the
% global cells number increase.
for i = 1 : 8
    boostcells = sum(sum((bone == site.tumor) & maschere_ott(:,:,i))) - sum(sum((BONE(:, :, hour - 5) == site.tumor) & maschere_ott(:, :, i)));  

    if boostcells >= 0
       Prob_ott(i)=boostcells / (sum(sum((bone == site.tumor) - (BONE(:, :, hour - 5) == site.tumor)))); 
    else
       Prob_ott(i) = 0;
    end
end  

%Random choice of the mask where to adress the new vessel.
rng('shuffle');
order = randperm(8); %Random permutation
Prob_ott_rand = Prob_ott(order); %Random permutation of the probility vector
maschere_ott_rand = maschere_ott(:, :, order); %Random permutation of the 8 masks
rng('shuffle');
dovecercare = zeros(rows, columns);
dovecercare = (bone == site.tumor & ~bw2); % sites that could be adressed for a new vessel (bw2 ? la maschera gialla quadrata che nasce dalla funzione creomaskvasi che fa in modo che non creascano nuovi vasi nella zona centrale del tumore che gi? ha una densit? maggiore)
[xxx,yyy] = find(dovecercare == 1); %possibili candidati ad essere centri del vaso

if size(xxx, 1) == 0 %If empty i search where there is simply tumor
   [xxx, yyy] = find(bone == site.tumor);
end 

nn = numel(xxx); %Number of elements
rng('shuffle');
iii = randperm(nn); 
xxx_rand = xxx(iii); %random permutation of possible centers 
yyy_rand = yyy(iii);

for i = 1:new_vess_number
    
    
%     Area_vessels   = sum(sum(vessels)) * site_dim * site_dim; %Area of the vessels
%     if(Area_vessels * 15 > curr_tumor_area) % Control: Atum must be > of 15 times Area of the vessels
%        break
%     end
    
    sum_prob = 0; % I sum the probabilities of each part till i reach the random number. Stochastic way to sellect which part will be targeted by angiogenesis 
    rng('shuffle');
    number = rand;
    for t = 1:8
        sum_prob = sum_prob + Prob_ott_rand(t);
        if sum_prob > number
           break
        end 
    end  
    MASCHERA_OTT = maschere_ott_rand(:, : ,t); %selected mask
    % I want to create a new vessel in the selected region and NOT outside the tumor mass  
    Flag_Vess_at_edge = 1; %Put it to high
    Flag_Hai_beccato_la_giusta_maschera = 1; %Put it to high
    tic %Control
    while (Flag_Vess_at_edge == 1 || Flag_Hai_beccato_la_giusta_maschera == 1)                
          rng('shuffle');
          randomNumber = randi(size(xxx, 1)); %I randomly select the center of the new vessel
          Cyy = xxx_rand(randomNumber); 
          Cxx = yyy_rand(randomNumber);
          a_vessel = random(pd_a_kernels); %Randomly select major semixis 
          b_vessel = random(pd_b_InverseGaussian); %Randomly select minor semixis 

          if b_vessel > a_vessel %if b>a I switch the 2 variables
             c_vessel = a_vessel;
             a_vessel = b_vessel;
             b_vessel = c_vessel;
          end 

          % Check on new center position
          if bone(Cyy + ceil(b_vessel / site_dim), Cxx) >= 2 && bone(Cyy - ceil(b_vessel / site_dim), Cxx) >= 2 && bone(Cyy, Cxx + ceil(a_vessel / site_dim)) >= 2 && bone(Cyy, Cxx - ceil(a_vessel / site_dim)) >= 2 %controlla i 4 vicini e vede se aggiungendo un possibile vaso si sfora la dimensione del tumore
             Flag_Vess_at_edge = 0;
          end 

          %If You have selected the corrected mask
          if MASCHERA_OTT(Cyy, Cxx)==1
             Flag_Hai_beccato_la_giusta_maschera = 0;
          end

          if  Flag_Vess_at_edge == 1 || Flag_Hai_beccato_la_giusta_maschera == 1
              clear  NewVesselCenter
          end

          tempotrascorso = toc;
          if tempotrascorso > 5  %If it's passed too much time not having found a solution i break the cycle
              break
          end 
    end  

    %New vessel vector: Cohordinates and its dimensions
    %if cabo therapy is in action i must add a column at
    %NewVesselsCenter to see whether the vessel is targeted
    %by cabo or not 

    if flag_cabo == 1
       NewVesselCenter = [Cyy, Cxx, a_vessel, b_vessel, 0]; %when cabo theraphy starts we add a column at center_vessels to see whether a vessel has been targeted or not by cabo. So we must add this column to New_center_vessels
       vessels_created_during_cabo = vessels_created_during_cabo + 1;
    elseif flag_cabo == 0
       NewVesselCenter = [Cyy, Cxx, a_vessel, b_vessel];
    end 

    %I update the vessels mask
    for jj = 1 : rows
        for  kk = 1 : columns
             if bone(jj, kk) == site.tumor  
                % Distance site-center of the vessel
                dist_centro_vaso = sqrt((X(jj,kk)-Cxx)^2 + (Y(jj,kk)-Cyy)^2);

                if dist_centro_vaso <= ellipse_radius(ceil(a_vessel / site_dim), ceil(b_vessel / site_dim), abs(X(jj, kk) - Cxx), abs(Y(jj, kk) - Cyy))

                   vessels(jj, kk) = 1; % Them the site will belong to the new vessel
                end
             end 
        end 
    end 
    
    Area_vessels   = sum(sum(vessels)) * site_dim * site_dim; %Area of the vessels
    center_vessels = [center_vessels; NewVesselCenter]; %I add the raw of the new vessel
    clear NewVesselCenter
    if(Area_vessels * 15 > curr_tumor_area) % Control: Atum must be > of 15 times Area of the vessels
       break
    end

    
end 
clear MASCHERA_OTT
%I update the bone matrix adding vessels
[xxxxx,yyyyy] = find((bone == site.tumor & vessels) == 1);

for kkk = 1 : size(xxxxx, 1)
    bone(xxxxx(kkk), yyyyy(kkk)) = site.vessel;
end

clear NewVesselCenter xxxxx  yyyyy
clear  xxx yyy iii dovecercare i xxx_rand  yyy_rand maschere_ott