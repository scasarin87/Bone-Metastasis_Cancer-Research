%% Thirty Hours Angiogenesis %%

%  Here we perform the new vessels generation every 20 hours

% Compute tumor area at the current iteration
curr_tumor_area = sum(sum(bone == site.tumor | bone == site.vessel | bone == site.vessel_cabo)) * site_dim * site_dim; %Area in squared µm           
% Compute the vessels to be generated with this angiogenetic step
new_vess_number = round((Mean_Vessel_Density) * curr_tumor_area) - vessels_number; 

% Vessels mask creation: finds the core of pre-established vessels; generates new ones outside of it
mask = ones(rows, columns); 
% Individuates a core of already established vessels. i wil try to generate new vessels outside it
vessels_mask_creation; 

rng('shuffle');
dovecercare = (bone == site.tumor & ~bw2);
[xxx,yyy] = find(dovecercare == 1);
if size(xxx, 1)==0
    [xxx,yyy] = find(bone == site.tumor);
end

nn=numel(xxx);
iii=randperm(nn);
xxx_rand=xxx(iii);
yyy_rand=yyy(iii);

for i = 1:new_vess_number
    
%     Area_vessels = sum(sum(vessels))*site_dim*site_dim; 
%     if(Area_vessels * 15 > curr_tumor_area)
%        break 
%     end 
    
    Flag_Vess_at_edge=1;
       while (Flag_Vess_at_edge==1)
            rng('shuffle');
            randomNumber = randi(size(xxx,1));
            Cyy  = xxx_rand(randomNumber);
            Cxx = yyy_rand(randomNumber);
            a_vessel=random(pd_a_kernels);
            b_vessel=random(pd_b_InverseGaussian);

%             if b_vessel > a_vessel
%                  c_vessel=a_vessel;
%                  a_vessel=b_vessel;
%                  b_vessel=c_vessel;
%             end

            if bone(Cyy+ceil(b_vessel/site_dim),Cxx)>=site.tumor && bone(Cyy-ceil(b_vessel/site_dim),Cxx)>=site.tumor && bone(Cyy,Cxx+ceil(a_vessel/site_dim))>=site.tumor && bone(Cyy,Cxx-ceil(a_vessel/site_dim))>=site.tumor
                Flag_Vess_at_edge=0;
            end
       end

       if flag_cabo == 1
          NewVesselCenter = [Cyy, Cxx, a_vessel, b_vessel,0];
          % vessels_created_during_cabo = vessels_created_during_cabo + 1;
       elseif flag_cabo == 0
              NewVesselCenter = [Cyy, Cxx, a_vessel, b_vessel];
       end

        for jj = 1 : rows
            for  kk = 1 : columns
               if bone(jj,kk) == site.tumor 
                   dist_centro_vaso = sqrt((X(jj,kk)-Cxx)^2 + (Y(jj,kk)-Cyy)^2);
                   if dist_centro_vaso <= ellipse_radius(ceil(a_vessel/site_dim),ceil(b_vessel/site_dim),abs(X(jj,kk)-Cxx),abs(Y(jj,kk)-Cyy))
                       vessels(jj, kk)=1;
                   end
               end
            end
        end
    Area_vessels = sum(sum(vessels))*site_dim*site_dim; 
    center_vessels = [center_vessels; NewVesselCenter];
    clear NewVesselCenter 
    if(Area_vessels * 15 > curr_tumor_area)
       break 
    end 
    
end

[xxxxx,yyyyy]=find((bone==site.tumor & vessels)==1);

for kkk=1:size(xxxxx,1)
    bone(xxxxx(i),yyyyy(i))=site.vessel;
end

clear NewVesselCenter xxxxx  yyyyy
clear  xxx yyy iii dovecercare i xxx_rand yyy_rand