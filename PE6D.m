% *******************************************************************************
% Rough pose estimation based on weak-perspective geometry model 
% Function: PE6D.m
% - Input: original image, pattern period T, 
%          camera nominal pixel size t & focallenght f
% - Output:alpha_hat,beta_hat,tz_hat
%

function [alpha_hat,beta_hat,tz_hat]=PE6D(data,T,t,f,window)
% initialization
    I=data; %input
    [N,M]=size(I);% number of pixels on image sensor 
    v0=ceil((N+1)/2);
    u0=ceil((M+1)/2); 
    p=f/t; % focal length (unit: pixel)
    const1=u0*t;%unit:pixel-->um
    const2=v0*t;%unit:pixel-->um
    
%% ------- Weak-perspective geometry model to estimate alpha,beta,tz ---------------

    vp=round(N/4);
    up=round(M/4);
    vC=vp; 
    uC=up; 
    vD=N-vp; 
    uD=up;
    vA=N-vp;
    uA=M-up;
    vB=vp; 
    uB=M-up;
    
 % extract phases of sub-images
    row_num=window;
    col_num=row_num; 
    row_start=vp-(row_num-1)/2;
    row_end=vp+(row_num-1)/2;
    col_start=up-(col_num-1)/2;
    col_end=up+(col_num-1)/2;

% calculate side length£ºCD,AB,BC,DA 
    [phase001,phase002]=PE2D_length1(N,I,row_num,row_start,row_end,col_start,col_end); 
    phase01=unwrap(phase001);
    phase02=unwrap(phase002);
    line01=(phase01(N-2*vp+1)-phase01(1))/(2*pi)*T;
    line02=(phase02(N-2*vp+1)-phase02(1))/(2*pi)*T;
    length_CD=sqrt(line01^2+line02^2); % CD

    [phase001,phase002]=PE2D_length2(N,I,row_num,row_start,row_end,col_start,col_end,M);
    phase01=unwrap(phase001);
    phase02=unwrap(phase002);
    line01=(phase01(N-2*vp+1)-phase01(1))/(2*pi)*T;
    line02=(phase02(N-2*vp+1)-phase02(1))/(2*pi)*T;
    length_AB=sqrt(line01^2+line02^2);  % AB

    [phase001,phase002]=PE2D_length3(M,I,col_num,row_start,row_end,col_start,col_end);
    phase01=unwrap(phase001);
    phase02=unwrap(phase002);
    line01=(phase01(M-2*up+1)-phase01(1))/(2*pi)*T;
    line02=(phase02(M-2*up+1)-phase02(1))/(2*pi)*T;
    length_BC=sqrt(line01^2+line02^2);  % BC

    [phase001,phase002]=PE2D_length4(M,I,col_num,row_start,row_end,col_start,col_end,N);
    phase01=unwrap(phase001);
    phase02=unwrap(phase002);
    line01=(phase01(M-2*up+1)-phase01(1))/(2*pi)*T;
    line02=(phase02(M-2*up+1)-phase02(1))/(2*pi)*T;
    length_DA=sqrt(line01^2+line02^2);  % DA

%% ------------------------------ Start points -------------------------

    tz1=(length_AB/const1+1)*f;
    tz2=(length_BC/const2+1)*f;  
    tz3=(length_CD/const1+1)*f;
    tz4=(length_DA/const2+1)*f;
    
    tz0=(tz1+tz2+tz3+tz4)/4; 
    alpha0=asin(2*(tz2-tz4)/(length_AB+length_CD));
    beta0=asin(2*(tz1-tz3)/(length_BC+length_DA));
    
%% -----------------------------Nonlinear optimization---------------------------------
  
  myfun=@(x)( (((uA-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (uB-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((vA-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (vB-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2-length_AB^2)^2+...
            (((uB-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (uC-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((vB-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (vC-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uB-u0)*sin(x(2)*beta0)-(vB-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2-length_BC^2)^2+...
            (((uC-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (uD-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((vC-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (vD-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uC-u0)*sin(x(2)*beta0)-(vC-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2-length_CD^2)^2+...;
            (((uD-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (uA-u0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((vD-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (vA-v0)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2+...
             ((-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uD-u0)*sin(x(2)*beta0)-(vD-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0))-...
              (-p)*(x(3)*tz0)*cos(x(1)*alpha0)*cos(x(2)*beta0)/((uA-u0)*sin(x(2)*beta0)-(vA-v0)*sin(x(1)*alpha0)*cos(x(2)*beta0)-p*cos(x(2)*beta0)*cos(x(1)*alpha0)))^2-length_DA^2)^2 );

    [x]=lsqnonlin(myfun,[1,1,1],[],[]);  

    % alpha_hat,beta_hat,tz_hat
    alpha_hat=x(1)*alpha0;  %unit:rad
    beta_hat=x(2)*beta0;    %unit:rad
    tz_hat=x(3)*tz0; %unit:um

end % end of PE6D.m
