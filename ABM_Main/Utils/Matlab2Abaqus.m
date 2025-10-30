%% Matlab2Abaqus %%

%  This function writes the .txt file with the cortical bone coordinates

function Matlab2Abaqus(nodes, filename)
    % File opening
    fileID = fopen(filename, 'w');
    %Generate nodes in input file
    [n_nodes, dim] = size(nodes);
    if dim == 2                              
        for node = 1:1:n_nodes
            fprintf(fileID, [num2str(nodes(node, 1)) ' ' num2str(nodes(node, 2)) ' ' '0' '\n']); 
        end 
    elseif dim == 3 
        for node = 1:1:n_nodes
            fprintf(fileID, [num2str(node) ', ' num2str(nodes(node,1)) ', ' num2str(nodes(node,2)) ', ' num2str(nodes(node,3)) '\n']); 
        end    
    end     
    fprintf(fileID, '\n');
    fclose(fileID);
end 