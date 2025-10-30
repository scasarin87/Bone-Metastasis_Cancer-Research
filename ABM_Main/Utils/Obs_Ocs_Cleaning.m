%This function is used to clean isolated obs/ocs
[row_obs,col_obs] = find(bone == site.osteoblast | bone == site.osteoclast);
for obs = 1:length(row_obs)

        
            if mod(col_obs(obs),2)==0
                j_k=2;
            else
                j_k=1;
            end
        

            liste=zeros(1,6);

            % look for neighbors and record their position
            for n3=1:6

                % (j3,k3) --> coordinates of the temporary (j,k) neighbor
                j3 = row_obs(obs)+directionx(n3,j_k);
                k3 = col_obs(obs)+directiony(n3);
                
                   

                % if obs/ocs are not attached with a CB site, then they are
                % replaced by tumor
                if bone(j3,k3) == site.tumor || bone(j3,k3) == site.tumor_edge || bone(j3,k3) == site.vessel || bone(j3,k3) == site.osteoblast || bone(j3,k3) == site.osteoclast
                    liste(1,n3) = 1;
                end
                
            end
            
                if liste(1,:) == 1
                   bone(row_obs(obs), col_obs(obs)) = site.tumor;
                end
            
end

% Update cortical bone_area
cortical_bone_sites(1, hour) = sum(sum(bone == site.cortical_bone)); 



 
