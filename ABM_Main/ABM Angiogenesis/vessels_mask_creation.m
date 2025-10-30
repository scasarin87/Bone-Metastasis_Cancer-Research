%voglio trovare il quadrato dove sta il tumore e traslarlo al centro per
%fare la maschera dei vasi;

%questo algoritmo funziona meglio se il tumore è centrale, quindi se il
%tumore è sul bordo prima lo trasla in centro, dopo di che crea la maschera
%e poi riporta il tumore in posizione iniziale

%I want to find the mask of the established vessels.
%I find the block of the only vessels;
% NN = round(10 + (5 / 12) * (rows - 30)); % Dimension of a box mask that will be then used to identify the vessels mask. 
NN = round(10 + (5 / 12) * (columns - 30));  % Its dimension its linearly dependent (Heuristically determined) to the
                                          % dimension Nx of thr grid.
[righe, colonne] = find(vessels == 1);
maschera_vasi = vessels(min(righe) : max(righe), min(colonne) : max(colonne));
stats = regionprops(vessels); % get the geometrical properties of the mask
r_c = fix(stats.Centroid(1)); % baricenter coordinates
c_c = fix(stats.Centroid(2));
maskvessels2 = zeros(rows, columns); %It will contains the vessel mask but traslated to the center of the matrix
% maskvessels2((round(central_col) - (c_c-min(righe))) : (round(central_row) + (max(righe) - c_c)), (round(central_col) - (r_c - min(colonne))) : (round(central_col) + (max(colonne) - r_c))) = maschera_vasi;
maskvessels2((round(central_row) - (c_c-min(righe))) : (round(central_row) + (max(righe) - c_c)), (round(central_col) - (r_c - min(colonne))) : (round(central_col) + (max(colonne) - r_c))) = maschera_vasi;
mask = zeros(rows, columns); %mask used to find the established vessels block
mask(NN : (rows - NN), NN : (columns - NN)) = 1; 
bw = zeros(rows, columns);
bw = activecontour(maskvessels2, mask); %getting the established vessel block Through this algorith that basically starting from the background it find the mask of the object in foreground
bw2 = zeros(rows, columns); %I re-traslate the mask bw in the original postion of the tumor
bw2(min(righe) : max(righe), min(colonne) : max(colonne)) = bw((round(central_row) - (c_c-min(righe))) : (round(central_row) + (max(righe) - c_c)), (round(central_col) - (r_c-min(colonne))) : (round(central_col) + (max(colonne) - r_c)));

clear righe colonne


