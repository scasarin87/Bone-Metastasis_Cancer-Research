% ellipse radius: given the length of the semi-axes of the original ellipse
% and the coordinate of point of the current site, the code gives in output
% the radius of the ellipse corresponding to the angle individuated by the
% couple of coordinates

function radius = ellipse_radius(a,b,X,Y) % a e b sono i semiassi maggiore e minore dell'ellisse, X e Y sono i valori del punto, da questi valori in ingresso ottengo il raggio
theta=atan(Y/X);
if X==0 && Y==0
    radius=0;
    theta=0;
else 
    radius = (a*b)/sqrt((b*cos(theta)).^2 + (a*sin(theta)).^2);
    theta=atan(Y/X);
end
end
