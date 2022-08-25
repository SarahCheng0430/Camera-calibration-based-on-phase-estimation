% *************************************************************************
% Phase Estimation method 
% - Function: PE2D.m
% - Usage: calculate the image center of 2D periodic pattern
% - Input: image, pattern pitch T,window function
% - Output:phaseX,phaseY,
% *************************************************************************

function [phase1, phase2, u1, v1, u2, v2, row1, col1, row2, col2, theta, Spectrum]=PE2D(data,Tx,Ty,win1,win2)
    data=data-mean(data(:)); 
    [N,M]=size(data);
    midpoint=ceil([(N+1)/2 (M+1)/2]); 
     if nargin<4
        win2=@hann;
        win1=@hann;
     elseif nargin==4
         win2=win1;
     end
    win=window2(N,M,win1,win2); % call window.m
    
    Spectrum=fftshift(fft2(data.*win)); %2D-FFT
    S_value0=abs(Spectrum); % Amplitude spectrum

    %% Find a pair of maximum amplitude points on S_value1
    S_value1=S_value0;  
    [row1,col1]=find(S_value1==max(max(S_value1)));
    u1=row1(2)-midpoint(1);
    v1=col1(2)-midpoint(2);
    %%%%%%%%%%%%%%%%%%%
    if row1(2)> midpoint(1) 
        point1=2;  
        k1=(u1/N)/(v1/M);
        theta1=atan(1/k1)/pi*180; % unit£ºdegree
    else  
        S_value1(row1(1)-2:row1(1)+2,col1(1)-2:col1(1)+2)=0;
        S_value1(row1(2)-2:row1(2)+2,col1(2)-2:col1(2)+2)=0;
        
        [row1,col1]=find(S_value1==max(max(S_value1)));
        u1=row1(2)-midpoint(1);
        v1=col1(2)-midpoint(2);  
        point1=2;  
        k1=(u1/N)/(v1/M);
        theta1=atan(1/k1)/pi*180; % unit£ºdegree
    end    
        
    %% Find the other pair of maximum amplitude points on S_value2
    S_value2=S_value0; 
    S_value2(row1(1)-2:row1(1)+2,col1(1)-2:col1(1)+2)=0;
    S_value2(row1(2)-2:row1(2)+2,col1(2)-2:col1(2)+2)=0;
    [row2,col2]=find(S_value2==max(max(S_value2)));
    u2=row2(2)-midpoint(1);
    v2=col2(2)-midpoint(2);
    point2=2; 
    if u2>0
        u2=row2(1)-midpoint(1);
        v2=col2(1)-midpoint(2);
        point2=1;
    end
    k2=(u2/N)/(v2/M);
    theta2=atan(1/k2)/pi*180;
    
    %
    while abs(theta1-theta2)<80 || abs(theta1-theta2)>100
        temp11=max(row2(1)-2,1);
        temp12=min(row2(1)+2,N);
        temp13=max(col2(1)-2,1);
        temp14=min(col2(1)+2,M);
        temp21=max(row2(2)-2,1);
        temp22=min(row2(2)+2,N);
        temp23=max(col2(2)-2,1);
        temp24=min(col2(2)+2,M);
        S_value2(temp11:temp12,temp13:temp14)=0;
        S_value2(temp21:temp22,temp23:temp24)=0;
        [row2,col2]=find(S_value2==max(max(S_value2)));
        u2=row2(2)-midpoint(1);
        v2=col2(2)-midpoint(2);
        point2=2; 
        if u2>0
            u2=row2(1)-midpoint(1);
            v2=col2(1)-midpoint(2);
            point2=1;  
        end
        k2=(u2/N)/(v2/M);
        theta2=atan(1/k2)/pi*180;  
    end

    %% Calculate PhaseX and PhaseY based on phase estimation method
    phase1=angle(Spectrum(row1(point1),col1(point1)))+pi*u1*(N-1)/N+pi*v1*(M-1)/M;
    phase1=mod(phase1,2*pi);
    phase2=angle(Spectrum(row2(point2),col2(point2)))+pi*u2*(N-1)/N+pi*v2*(M-1)/M;
    phase2=mod(phase2,2*pi);
    
    %% theta
    if theta1<0
        theta=180+theta1;
    else
        theta=theta1;
    end   
    
end
