function [dclb] = Dclb(clb,G,Fy,m,v)
%������΢��
dclb=(G*cos(clb*pi/180)-Fy)/(m*v);
end

