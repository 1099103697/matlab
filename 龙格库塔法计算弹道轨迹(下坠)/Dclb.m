function [dclb] = Dclb(clb,G,Fy,m,v)
%ÅÀÉý½ÇÎ¢·Ö
dclb=(G*cos(clb*pi/180)-Fy)/(m*v);
end

