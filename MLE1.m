% ***********************************************************************************
% - Iteratively optimize the calibration parametes 
% - Cost function based on reprojection errors
% ***********************************************************************************
% REF:	   "A Flexible New Technique for Camera Calibration"
%           - Zhengyou Zhang 
%           - Microsoft Research 

% cost function MLE1

function f = MLE1(H, m, M)
% unpack the params
    
    %h=[H([1:3]); H([4:6]); H([7:8]),1];
    H=reshape(H,3,3);
    h=H;
% unpack the params_const
    
    m=[m([1:2],:); ones(1,size(m,2))];
    M=[M([1:2],:); ones(1,size(M,2))];
    
    X=h*M;
    X=[X(1,:)./X(3,:) ; X(2,:)./X(3,:); X(3,:)./X(3,:)];
    
    res=m-X;
    req=[res(1,:), res(2,:)]; 
    
    
    f = req;
    
end % end of function