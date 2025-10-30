% MASCHERE IN 8 Quadranti;
function maschere_ott=octants_mask_creation(r_c, c_c, rows, columns)
    
    % FARE UN CICLO FOR AL POSTO DI STO SCHIFO
    
    Cy = r_c; %Coordinates of the center of the 8 masks
    Cx = c_c;
    I1 = zeros(rows, columns); %Equation of the 1° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;
             if  y < x && x > 0 && y > 0
                 I1(i, j) = 1;
             end
        end
    end
    
    I2 = zeros(rows, columns); %Equation of the 2° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;
             if  y > x && x > 0 && y > 0
                 I2(i, j) = 1;
             end
        end
    end
    
    I3 = zeros(rows, columns); %Equation of the 3° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
            y = Cy - i;
            x = j - Cx;
            if  y < -x && x < 0 && y > 0
                I3(i,j)=1;
            end
        end
    end
    
    
    I4 = zeros(rows, columns); %Equation of the 4° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;
             if  y > -x && x < 0 && y > 0
                 I4(i, j) = 1;
             end
        end
    end
    
    I5 = zeros(rows, columns); %Equation of the 5° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;
             if  y > x && x < 0 && y < 0
                 I5(i, j) = 1;
             end
        end
    end
    
    I6 = zeros(rows, columns); %Equation of the 6° octant based on the bisector equation
    for i = 1 : rows
        for j = 1 : columns
            y = Cy - i;
            x = j - Cx;
            if  y < x && x < 0 && y < 0
                 I6(i, j) = 1;
            end 
        end
    end
    
    I7 = zeros(rows,columns); %Equation of the 7° octant based on the bisector equation
    for i=1:rows

        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;

             if  y < -x && x > 0 && y < 0
                 I7(i, j) = 1;
             end
        end
    end
    % imagesc(I7)

    I8 = zeros(rows, columns); %Equation of the 8° octant based on the bisector equation
    for i = 1 : rows

        for j = 1 : columns
             y = Cy - i;
             x = j - Cx;

             if  y > -x && x > 0 && y < 0
                 I8(i, j)=1;
             end
        end
    end
    % imagesc(I8)
    maschere_ott=zeros(rows,columns,8);
    maschere_ott(:,:,1)=I1;
    maschere_ott(:,:,2)=I2;
    maschere_ott(:,:,3)=I4;
    maschere_ott(:,:,4)=I3;
    maschere_ott(:,:,5)=I5;
    maschere_ott(:,:,6)=I6;
    maschere_ott(:,:,7)=I7;
    maschere_ott(:,:,8)=I8;
    clear I1 I2 I3 I4 I5 I6 I7 I8
end