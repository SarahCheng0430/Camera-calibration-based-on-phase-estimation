% *************************  Remapping  *********************************
% Function :remapping.m
%   - Usage:perspective correction, remapping to a fronto-parallel plane
% - Input: original image; tilt angle estimation parameters
% - Output:new fronto-parallel image 

function [imgnew,delta_u,delta_v,H]=remapping(img,alpha_degree,beta_degree,In)
[M,N] = size(img);
theta_degree=0;%set theta==0
alpha=alpha_degree/180*pi;% unit:rad
beta=beta_degree/180*pi;
theta=theta_degree/180*pi; 

% R=Rx*Ry*Rz
Rx=[1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Ry=[cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
R=Rz*Ry*Rx; 

%homography;
H=In*inv(R)*inv(In); 
% a11=H(1,1);a12=H(1,2);a13=H(1,3);
% a21=H(2,1);a22=H(2,2);a23=H(2,3);
a31=H(3,1);a32=H(3,2);a33=H(3,3);

%% Construct a new image
% corners for new image
Pixel1=H*[1 1 1]'/(a31*1+a32*1+a33);  
Pixel2=H*[1 N 1]'/(a31*1+a32*N+a33);  
Pixel3=H*[M 1 1]'/(a31*M+a32*1+a33);  
Pixel4=H*[M N 1]'/(a31*M+a32*N+a33);  

height=round(max([Pixel1(1) Pixel2(1) Pixel3(1) Pixel4(1)])-min([Pixel1(1) Pixel2(1) Pixel3(1) Pixel4(1)]));     
width=round(max([Pixel1(2) Pixel2(2) Pixel3(2) Pixel4(2)])-min([Pixel1(2) Pixel2(2) Pixel3(2) Pixel4(2)]));      
% disp('imgremapped size');
% disp([height, width]);

imgnew=zeros(height,width);
delta_u=round((min([Pixel1(1) Pixel2(1) Pixel3(1) Pixel4(1)])));   %Offset u 
delta_v=round((min([Pixel1(2) Pixel2(2) Pixel3(2) Pixel4(2)])));   %Offset v


% (u2,v2)-(height width)-----(u1,v1)-(M N)
for u2 = 1+delta_u:height+delta_u       
    for v2 = 1+delta_v:width+delta_v
       
        pix0=inv(H)*[u2 v2 1]';  
        pix=pix0/pix0(3);        

        if pix(1)>=1 && pix(2)>=1 && pix(1)<=M && pix(2)<=N
            %Bilinear interpolation
            u1=pix(1);v1=pix(2);
            u11=floor(u1);v11=floor(v1);u12=u11+1;v12=v11+1;
            imgnew(u2-delta_u,v2-delta_v)=(u12-u1)*(v12-v1)*img(u11,v11)+(u12-u1)*(v1-v11)*img(u11,v12)+(u1-u11)*(v12-v1)*img(u12,v11)+(u1-u11)*(v1-v11)*img(u12,v12);    
        end  
    end
end


end



