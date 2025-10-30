%% Processing Bone Image
%  This function is used to perform bone image preprocessing according to
%  the following steps:
%  Binarization: to get only 1 and 0s from the original image
%  Reshaping   : to fit in vivo real dimensions (x_real_dim to be defined
%                acoordingly, in micrometers)
%  Padding     : I'll perform padding either horizontally or vertically to 
%                have a final square shape 
%  Finally, we define the hexagonal grid.

% Input : image           -> the bone image which has just been loaded
%         x_real_dim      -> real dimension on x_axis[in um] (ask Eleonora)
%         site_dim        -> ABM rescale factor, currently set at 20.833um
%
% Output: processed_image -> the image is processed (binarization, reshape, padding)
%         rows, columns   -> size of the processed image
%         X, Y            -> hexagonal grid created according to the matrix dimensions

function [processed_image, rows, columns, X, Y] = processing_BoneImage(image, x_real_dim, site_dim)
    %% Binarization: mask creation
    binary_image = imbinarize(image, 'global'); % Binarization
    binary_image = binary_image(:, :, 1);       % I just need 2D binarized image
    
    %% Reshape: mask must fit in vivo dimensions
    curr_x_dim = size(binary_image, 2);    % Get current x image shape
    curr_x_dim_um = curr_x_dim * site_dim; % Get the x shape in um
    
    if curr_x_dim_um > x_real_dim % I need the reshape factor to be < 1
        reshape_factor = x_real_dim / curr_x_dim_um;
    else
        reshape_factor = curr_x_dim_um / x_real_dim;
    end
    
    % Image reshape to fit in vivo bone dimensions
    reshaped_image = imresize(binary_image, reshape_factor); 
    
    %% Padding: I want x and y dimensions to be the same
    [new_y_dim, new_x_dim] = size(reshaped_image);
    
    % Padding is performed by adding a number of rows or colums to fill the
    % gap between the x and y size, obtaining a final image shape where
    % x_dim = y_dim
        
    if new_x_dim >= new_y_dim
        % I'll pad horizzontally
        diff = new_x_dim - new_y_dim;
        pad_upvector_dim     = floor(diff / 2); % floor and ceil are used to deal with numbers which can't be divided by 2
        pad_downvector_dim   = ceil(diff / 2);
        up_pad_vector        = ones(pad_upvector_dim, new_x_dim);
        down_pad_vector      = ones(pad_downvector_dim, new_x_dim);
        pad_image = [up_pad_vector; reshaped_image; down_pad_vector];    
        
    elseif new_y_dim > new_x_dim
        % I'll pad vertically
        diff = new_y_dim - new_x_dim;
        pad_leftvector_dim   = floor(diff / 2); % floor and ceil are used to deal with numbers which can't be divided by 2
        pad_rightvector_dim  = ceil(diff / 2);
        left_pad_vector      = ones(new_y_dim, pad_leftvector_dim);
        right_pad_vector     = ones(new_y_dim, pad_rightvector_dim);
        pad_image = [left_pad_vector reshaped_image right_pad_vector];            
    end
    
    processed_image = pad_image; 
    
    clear pad_image reshaped_image bone_image
    
    %% Hexagonal Grid definition 
    %  I create the hexagonal grid according to the size of the new image
    [rows, columns] = size(processed_image);
    X = zeros(rows, columns); Y = zeros(rows, columns); %Y refers to rows, X to columns 
    
    for jj=1:rows
        for kk=1:columns

            if mod(kk,2)==0
                Y(jj,kk)=2*jj; X(jj,kk)=2*kk;
            else
                Y(jj,kk)=2*jj+1; X(jj,kk)=2*kk;
            end

        end
    end
    
    clear jj kk
    X=0.5.*X; Y=0.5.*Y;
    
end
