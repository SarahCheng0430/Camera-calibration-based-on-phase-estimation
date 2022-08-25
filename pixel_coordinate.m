% ************************* pixel coordinate *********************************
% - Usage:Pixel reversion transform of Control points
% - Input: [u2,v2] in cropped remapped image, H: homography
% - Output:[u1,v1] in original image

function [u1,v1]=pixel_coordinate(u2,v2,H,delta_u,delta_v)
%(u2,v2)-->(u1,v1)
u2=u2+delta_u;
v2=v2+delta_v;
%Pixel reversion
pix0=inv(H)*[u2 v2 1]'; %H homography:x'=H*x
pix=pix0/pix0(3);

u1=pix(1);
v1=pix(2);

end