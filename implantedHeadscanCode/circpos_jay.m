function [angle] = circpos_jay(posx, posy)

[xc,yc,~,~] = circfit(posx,posy);

ax = posx-xc;
ay = posy-yc;

angle = NaN(1,size(posx,2));
for aa = 1:size(posx,2)
    x = ax(aa);
    y = ay(aa);
    
    %first quadrant
    if x>0 && y>=0 
        angle(aa) = asind(y/sqrt(x^2+y^2));
    %second quadrant    
    elseif x<=0 && y>0
        angle(aa) = 180-asind(y/sqrt(x^2+y^2));
    %third quadrant
    elseif x<0 && y<=0
        angle(aa) = asind(abs(y)/sqrt(x^2+y^2))+180;
    %fourth quadrant
    elseif x>=0 && y<0
        angle(aa) = 360-asind(abs(y)/sqrt(x^2+y^2));
    end
    
end
angle = angle+180;
angle(angle>360)=angle(angle>360)-360;
