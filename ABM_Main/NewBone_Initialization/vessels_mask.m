%% Vessels Mask

% Input  -> rows, columns   : are x and y size of the processed image matrix
%           tumor           : is the tumor mask (= 1 if site is PCa cell)
%           X, Y            : hexagonal grid sites
%           site_dim        : pixel -> um conversion rate (set to 500/24)
%
% Output -> vessels         : 1 -> vessel site
%                             0 -> other                                                                 
%           center_vessels_0: vessels' center matrix with axis dimension

function [vessels, center_vessels_0] = vessels_mask(rows, columns, tumor, site_dim, X, Y)

    load VesselsProperties
    
    %I compute the number of vessels that have to be created and I randomly select their size and centers  
    Area_tumor        = sum(sum(tumor == 1)) * site_dim * site_dim; %Area in squared um      
    Number_of_vessels = ceil((Mean_Vessel_Density) * Area_tumor); %Number of initial vessels to be created

    % Poisson disc sampling: The distance between vessels must be at least of a number of pixel determined by spacing variable
    sizeI   = [rows, columns]; 
    vessels = zeros(rows, columns);
    spacing = 3;  
    pts     = poissonDisc(sizeI, spacing); % It Determines the poissons points (sampling to define vessels position)
    II      = zeros(rows, columns);        %Mask of the poisson points

    for i = 1:size(pts,1)
        II(round(pts(i,1)), round(pts(i,2))) = 1;
    end
    
    clear i

    %I will draw the centers of the vessels among the points belonging to tumor and being poisson points.
    %I perform a random permutation of the points list to determine the centers of the vessels.
    [xx, yy] = find(tumor == 1 & II); 
    rng('shuffle')
    n       = numel(xx); 
    ii      = randperm(n);
    xx_rand = xx(ii);
    yy_rand = yy(ii);
    possible_centers = [xx_rand,yy_rand]; %List of possibile center coordinates
    center_vessels   = []; % Each row corresponds to a vessel
                           % 1° and 2° columns: cohordinate of the vessel center 
                           % 3° and 4° columns: Semiaxis of the vessel dimension 

    %Cycle that will create the vessels
    for iiii = 1:Number_of_vessels                

        %Defining randomly vessels center and semiaxis
        %I select the first number of vessels center from NewVesselCenter Flag_Vess_at_edge = 1;

        %Randomly select new vessel center
        Cyy  = possible_centers(iiii,1); 
        Cxx = possible_centers(iiii,2);

        %Randomly select new vessel semi-axis
        a_vessel=random(pd_a_kernels);
        b_vessel=random(pd_b_InverseGaussian);

        %a_vessel always > b_vessel
        if b_vessel > a_vessel
           c_vessel=a_vessel;
           a_vessel=b_vessel;
           b_vessel=c_vessel;
        end 

        NewVesselCenter = [Cyy, Cxx, a_vessel, b_vessel];
        center_vessels = [center_vessels; NewVesselCenter]; %center vessels matrix update
        clear NewVesselCenter  

        %Vessel Mask will be updated with respect to the vessels that are being created                            
        %The vessels mask will be = 1 where a new vessel has been created (0 otherwise)

        for jj = 1:rows
             for  kk = 1:columns
                  if tumor(jj,kk) == 1 
                     dist_centro_vaso = sqrt((X(jj,kk)-Cxx)^2 + (Y(jj,kk)-Cyy)^2); %I compute the site-vessel center distance
                     if dist_centro_vaso <= ellipse_radius(ceil(a_vessel/site_dim),ceil(b_vessel/site_dim),abs(X(jj,kk)-Cxx),abs(Y(jj,kk)-Cyy))
                        vessels(jj,kk) = 1; %Vessels mask
                     end  
                  end 
             end 
        end
         
        clear jj kk

        vessels(Cyy, Cxx) = 1; %minimum dimension of a vessel is 1 pixel 

        %Control of a possibile Over-creation of vessels 
        Area_vessels = sum(sum(vessels))*site_dim*site_dim; % in µm^2        
        if(Area_vessels*15 > Area_tumor) %Parameter of control: Atum>Avess*15 *from experimental images*
           Number_of_vessels = iiii; %I stop iterations if the area criterion is not met
           break 
        end                

    end
    
    clear iiii
    
    center_vessels_0 = center_vessels; %Starting vessels matrix
    
end