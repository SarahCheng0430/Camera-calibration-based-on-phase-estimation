function [CameraParams] = main_calibration_app(imgs, option, opt_value, T, t, f)
    
    command = "";
    if option == 0
        command = "normal";
    elseif option == 1
        command = "blur with r = " + num2str(opt_value(1));
    elseif option == 2
        command = "noise with [mu, sigma] = [" + num2str(opt_value(1)) + " , " + num2str(opt_value(2)) + "]";
    end

warning('off');
%% initialize parameters
    [M, N, num] = size(imgs);
    u0=((M+1)/2);
    v0=((N+1)/2);
    ctr=[u0,v0]; % original image center
    win_size=141;
    dX=t;dY=t; % pixel size
    fx=f/dX; %normalized focal length(unit:pixel)
    fy=f/dY; 
    skew=0; %skew between u/v axis
    In=[fx skew u0;%initialized Intrinsic matrix
        0 fy v0;
        0 0 1];

%% Phase estimation method to calculate point correspondence
    count=0; % count the num of valid images
    
    for j =1:1:num   % input images 
        img = process_img(imgs(:,:,j), option, opt_value);
        imgs(:,:,j) = img;
        
        % step 1 preprocess:original image-->cropped imgnew
        [imgnew,valid] = preprocess(img,ctr,win_size);
        
    % valid is used to judge this image is valid or not
        if valid  
            count=count+1;
            
            % step 2 estimate the tilt angle using cropped imgnew
            [beta0,alpha0,tz0]=PE6D(imgnew,T,t,f,win_size);
            alpha_degree=alpha0*180/pi; %unit:degree
            beta_degree=beta0*180/pi; %unit:degree
            tz_um=tz0; %unit:um

            % step 3 remapping original image --> fronto-parallel image
            [imgremapped,delta_u,delta_v,H]=remapping(img,alpha_degree,beta_degree,In);
%             figure
%             imshow(imgremapped);

            % step 4  find the pattern center marker Ow in the remapped image
            % [uOw,vOw] is the image coordinate of pattern origin Ow in the remapped image
            [uOw,vOw,imgcut] = find_center(imgremapped);

            % step 5 calculate the 4 point pixel coordinate in original image and its world coordinate
            [Q,P]=cali_coordinate(imgcut,uOw,vOw,delta_u,delta_v,H,T);
            q(:,:,count)=Q;%pixel plane points 2*npts*num (size=2*8*valid_num)
            p(:,:,count)=P;%object plane points 2*npts*num (size=2*8*valid_num)
            
            disp('done');
            clc;
            close all;
            
        else 
            % if valid=false£¬then next image
            continue
        
        end
                 
    end

%% Calculate intrinsic and distortion parameters
% 1. Transfer to homogeneous coordinates ---------------------------------
    [~,npts,num]=size(p);    %rows=2 refers to [X,Y] world coordinate
    matrixone=ones(1,npts);

    for i=1:num
        q(3,:,i)=matrixone; %m=[u,v]' ---> m=[u,v,1]' 
    end
    for i=1:num
        p(3,:,i)=matrixone; % M=[X,Y]' ---> M=[X,Y,1]' 
    end

% 2. Estimate the H ------------------------------------------
% H=A*[r1,r2,t]
    for i=1:num %transfer M->M(:,:,i)
        H(:,:,i)=homography2d(p(:,:,i),q(:,:,i))';%%%homography2d.m 
    end
    
% 3. Calculate intrinsic parameters --------------------------------
    V=[];
    for flag=1:num %Eq(7)
        v12(:,:,flag)=[H(1,1,flag)*H(2,1,flag), H(1,1,flag)*H(2,2,flag)+H(1,2,flag)*H(2,1,flag), H(1,2,flag)*H(2,2,flag), H(1,3,flag)*H(2,1,flag)+H(1,1,flag)*H(2,3,flag), H(1,3,flag)*H(2,2,flag)+H(1,2,flag)*H(2,3,flag), H(1,3,flag)*H(2,3,flag)];
        v11(:,:,flag)=[H(1,1,flag)*H(1,1,flag), H(1,1,flag)*H(1,2,flag)+H(1,2,flag)*H(1,1,flag), H(1,2,flag)*H(1,2,flag), H(1,3,flag)*H(1,1,flag)+H(1,1,flag)*H(1,3,flag), H(1,3,flag)*H(1,2,flag)+H(1,2,flag)*H(1,3,flag), H(1,3,flag)*H(1,3,flag)];
        v22(:,:,flag)=[H(2,1,flag)*H(2,1,flag), H(2,1,flag)*H(2,2,flag)+H(2,2,flag)*H(2,1,flag), H(2,2,flag)*H(2,2,flag), H(2,3,flag)*H(2,1,flag)+H(2,1,flag)*H(2,3,flag), H(2,3,flag)*H(2,2,flag)+H(2,2,flag)*H(2,3,flag), H(2,3,flag)*H(2,3,flag)];
        V=[V;v12(:,:,flag);v11(:,:,flag)-v22(:,:,flag)]; %Eq(8)
    end
    
    k=V'*V;        %b=0
    [~,~,d]=svd(k);%[u,s,v]=svd(A),A=USV'
    [~,~]=eig(k);  %Eigenvector [V,D]=eig(A)
    b=d(:,6);      
    
    v0=(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2);
    s=b(6)-(b(4)^2+v0*(b(2)*b(4)-b(1)*b(5)))/b(1);
    alpha_u=sqrt(s/b(1));
    alpha_v=sqrt(s*b(1)/(b(1)*b(3)-b(2)^2));
    skewness=-b(2)*alpha_u*alpha_u*alpha_v/s;
    u0=skewness*v0/alpha_u-b(4)*alpha_u*alpha_u/s;
    
    % intrinsic matrix 
    A=[alpha_u skewness u0
        0      alpha_v  v0
        0      0        1];
    
% 4.Calculate extrinsic matrix and distortion  -------------------------------------------
    D=[];
    d=[];
    Rm=[];
    for flag=1:num
        s=(1/norm(inv(A)*H(1,:,flag)')+1/norm(inv(A)*H(2,:,flag)'))/2;
        rl1=s*inv(A)*H(1,:,flag)';%r1=s*inv(A)*h1
        rl2=s*inv(A)*H(2,:,flag)';%r2=s*inv(A)*h2
        rl3=cross(rl1,rl2);       %r3=r1Xr2; % C = cross(A,B) returns the cross product of the vectors A and B. 
        % That is, C = A x B.  A and B must be 3 element vectors.
        RL=[rl1,rl2,rl3];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [U,~,V] = svd(RL);
        RL=U*V';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TL=s*inv(A)*H(3,:,flag)';%TL is Translation vector t
        RT=[rl1,rl2,TL];%H=A[r1 r2 t]
        XY=RT*p(:,:,flag); % p is points on model plane and XY refers to camera coordinate
        UV=A*XY;      %sm=A[R t]M,UV refers to pixel coordinate
        UV=[UV(1,:)./UV(3,:); UV(2,:)./UV(3,:); UV(3,:)./UV(3,:)]; % [u,v,1]
        XY=[XY(1,:)./XY(3,:); XY(2,:)./XY(3,:); XY(3,:)./XY(3,:)]; % [x,y,1]
        for j=1:npts
            D=[D; ((UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 )) , ((UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2) ; ((UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 )) , ((UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2) ];
            d=[d; (q(1,j,flag)-UV(1,j)) ; (q(2,j,flag)-UV(2,j))];
        end
        r13=RL(1,3);
        r12=RL(1,2);
        r23=RL(2,3);
        % Q1,Q2,Q3 are Euler angles
        Q1=-asin(r13);
        Q2=asin(r12/cos(Q1));
        Q3=asin(r23/cos(Q1));
        R_new=[Q1,Q2,Q3,TL'];
        Rm=[Rm, R_new];
    end
    
% Estimating radial distortion by alternation using function 
    k=inv(D'*D)*D'*d;
    
% Complete overall iterative refinement, using function 
    para=[Rm,k(1),k(2),alpha_u,skewness,u0,alpha_v,v0];
    % optimset Create/alter OPTIM OPTIONS structure.
    %options = optimset('LargeScale','off','LevenbergMarquardt','on'); 
    options = optimset('Display','off','TolX',eps,'TolFun',eps,'LargeScale','on',...
   'Algorithm','Levenberg-Marquardt');
    %lsqnonlin :Solves non-linear least squares problems¢˜
    %%%cost function MLE2
    [x,~,~,~,~]  = lsqnonlin( @MLE2, para, [],[],options, q, p);

%% output the result
    k = [x(num*6+1) x(num*6+2)];
    In = [x(num*6+3),x(num*6+4),x(num*6+5); 0,x(num*6+6),x(num*6+7); 0,0,1];
    Q = q(1:2,:,:);
    P = p(1:2,:,:);
    for ipa = 0:num-1
    	R_new = x([ipa*6+1:ipa*6+6]);
        Q1=R_new(1);
        Q2=R_new(2);
        Q3=R_new(3);
        TL=R_new([4:6])';
        RL=[cos(Q2)*cos(Q1)   sin(Q2)*cos(Q1)   -sin(Q1) ; 
            -sin(Q2)*cos(Q3)+cos(Q2)*sin(Q1)*sin(Q3)    cos(Q2)*cos(Q3)+sin(Q2)*sin(Q1)*sin(Q3)  cos(Q1)*sin(Q3) ;
            sin(Q2)*sin(Q3)+cos(Q2)*sin(Q1)*cos(Q3)    -cos(Q2)*sin(Q3)+sin(Q2)*sin(Q1)*cos(Q3)  cos(Q1)*cos(Q3)];
        RT=[RL(:,1:2) , TL];
        Rms(:,:,ipa+1) = RT;
        TmL(ipa+1,:) = rotationMatrixToVector(RL);
        Tm(ipa+1,:) = TL;
    end   
    
%     cam = cameraParameters('IntrinsicMatrix',In,'radialDistortion',[k(1) k(2)],'ImageSize',[M N], ...
%         'TranslationVectors',Tm,'RotationVectors',TmL,'WorldPoints',reshape(P(:,:,1),4,2));       
%     figure;showExtrinsics(cam, 'patterncentric');hold on;set(gca,'DataAspectRatio',[1 1 3]);    
%     figure;showExtrinsics(cam, 'cameracentric');hold on;set(gca,'DataAspectRatio',[1 1 3]);
    
% Calibration parameters output includes:
    CameraParams = struct('k', k, 'In', In, 'Rm', Rms, 'p', P, 'q', Q, 'opt', command);

    % process 3 types of images
    function imgn = process_img(image, opt, coef)
        if opt == 0
            imgn = mat2gray(image);
        elseif opt == 1
            r = coef(1);
            I = mat2gray(image);
            PSF = fspecial('disk',r);   %defocusing function
            imgn = imfilter(I,PSF,'symmetric','conv');  
        elseif opt == 2
            mu = coef(1); sigma = coef(2);
            I = image;
 %           imshow(I);
 %           imgn = mat2gray(imgn);
            [szM,szN] = size(I);
            R = normrnd(mu,sigma,szM,szN);
            image_noise = double(I) + R;
            imgn = uint8(round(image_noise));
            imgn = mat2gray(imgn);
        else
            disp('error option');
            imgn = 0;
        end
    end
end

