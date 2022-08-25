% ************************ cali_coordinate  *********************************
% Function cali_coordinate.m 
%   - Usage:obtain point correspondences based on phase estimation
%   - Call remapping.m 
%   - Call pixel_coordinate.m & world_coordinate.m
%   - Input: img, uO2,vO2,delta_u,delta_v,H,T
%   - Output:Qc,Pw

function [Qi,Pw]=cali_coordinate(imgcut,uO2,vO2,delta_u,delta_v,H,T)
%% -------------------Pixel coordinate (u,v,1)-----------------------------
% Calculate the pixel coordinates of the selected control points in the original image[u1,v1,1]
[M1,N1]=size(imgcut);
uO3=((M1+1)/2);
vO3=((N1+1)/2);
% O3=[uO3,vO3];

% [u3,v3] is coordinates in cropped remapped image
up=M1/4;vp=N1/4;
uA3=up;vA3=vO3; 
uB3=uO3;vB3=vp; 
uC3=M1-up;vC3=vO3;
uD3=uO3;vD3=N1-vp;
uE3=up;vE3=vp;
uF3=up;vF3=N1-vp;
uG3=M1-up;vG3=N1-vp;
uH3=M1-up;vH3=vp;

% (u2,v2)is coordinates in remapped image
shift_u=uO2-uO3;shift_v=vO2-vO3;
uA2=uA3+shift_u;vA2=vA3+shift_v;
uB2=uB3+shift_u;vB2=vB3+shift_v;
uC2=uC3+shift_u;vC2=vC3+shift_v;
uD2=uD3+shift_u;vD2=vD3+shift_v;
uE2=uE3+shift_u;vE2=vE3+shift_v;
uF2=uF3+shift_u;vF2=vF3+shift_v;
uG2=uG3+shift_u;vG2=vG3+shift_v;
uH2=uH3+shift_u;vH2=vH3+shift_v;

% call pixel_coordinate.m to obtain [u1,v1] 
[uA1,vA1]=pixel_coordinate(uA2,vA2,H,delta_u,delta_v);
[uB1,vB1]=pixel_coordinate(uB2,vB2,H,delta_u,delta_v);
[uC1,vC1]=pixel_coordinate(uC2,vC2,H,delta_u,delta_v);
[uD1,vD1]=pixel_coordinate(uD2,vD2,H,delta_u,delta_v);
[uE1,vE1]=pixel_coordinate(uE2,vE2,H,delta_u,delta_v);
[uF1,vF1]=pixel_coordinate(uF2,vF2,H,delta_u,delta_v);
[uG1,vG1]=pixel_coordinate(uG2,vG2,H,delta_u,delta_v);
[uH1,vH1]=pixel_coordinate(uH2,vH2,H,delta_u,delta_v);

% Pixel coordinates in input images
Pix_A=[uA1,vA1]';
Pix_B=[uB1,vB1]';
Pix_C=[uC1,vC1]';
Pix_D=[uD1,vD1]';
Pix_E=[uE1,vE1]';
Pix_F=[uF1,vF1]';
Pix_G=[uG1,vG1]';
Pix_H=[uH1,vH1]';
Qi=[Pix_A Pix_B Pix_C Pix_D Pix_E Pix_F Pix_G Pix_H]; %2*n 

%% ------------------ World coordinate (Xw,Yw,Zw=0,1)-------------------------
% Calculate world coordinate [Xw,Yw,Zw=0]   
[Aw,Bw,Cw,Dw,Ew,Fw,Gw,Hw]=world_coordinate(imgcut,T);
Pw=[Aw,Bw,Cw,Dw,Ew,Fw,Gw,Hw]; 

end %end of coordinate.m

