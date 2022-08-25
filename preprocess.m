function [imgnew,valid] = preprocess(img,ctr,win_size)
% ****************************************************************
% Function :preprocess.m
% - imgnew is used to calculate tilt angles in PE6D.m
% - judge this image is valid or not
% *****************************************************************

%Identify the area that contains the pattern in the img
% edge detection 
I = edge(img, 'log');
kernel = ones(10,10);
I_s = imdilate(I, kernel);
I_s = imerode(I_s, kernel);
I_s = imerode(I_s, kernel);
I_s = imdilate(I_s, kernel);

imLabel = bwlabel(I_s);
stats = regionprops(imLabel,'Area');    
area = cat(1,stats.Area);
index = find(area == max(area));        
imgnew = ismember(imLabel,index); 
% figure
% imshow(imgnew);title('Pattern Region')

% Extract 4 corner points containing pattern region
[h,w] = size(imgnew);
m1 = h+w;
m2 = 0;
n1 = h-1;
n2 = 1-w;
keypts = zeros(4,2);
for x = 1:h
    for y = 1:w
        if imgnew(x,y) == true
            sum_m = x + y;
            sum_n = x - y;
            if sum_m < m1
                keypts(1,:) = [x,y];
                m1 = sum_m;
            end
            if sum_m > m2
                keypts(2,:) = [x,y];
                m2 = sum_m;
            end
            if sum_n < n1
                keypts(3,:) = [x,y];
                n1 = sum_n;
            end
            if sum_n > n2
                keypts(4,:) = [x,y];
                n2 = sum_n;
            end
        end
    end
end
% disp(keypts);

%% crop a pattern area to estimate tilt angles
%ctr = [(h+1)/2, (w+1)/2]; 
%Symmetrically cut with the center point of the original image
function cut = find_cut(mainpts, center)
    function jud = judge(mainpts, pt)
        mainpts(:,3) = 0;
        pt(1,3) = 0;
        D = mainpts(1,:);
        B = mainpts(2,:);
        C = mainpts(3,:);
        A = mainpts(4,:);
        a = dot(cross(B-A, pt-A), cross(D-C,pt-C)) >= 0;
        b = dot(cross(A-D, pt-D), cross(C-B,pt-B)) >= 0;
        jud = a && b;
    end
    
    cut = 0;
    while cut <= 400      % set cropped image sidelength
        l = (cut+19)/2;
        j1 = judge(mainpts, center+[l,l]);
        j2 = judge(mainpts, center+[-l,l]);
        j3 = judge(mainpts, center+[l,-l]);
        j4 = judge(mainpts, center+[-l,-l]);
        if j1 && j2 && j3 && j4
            cut = cut + 20;  % 20 pixel step by step 
        else
            break;
        end
    end
    cut = (cut-1)/2;
end

cut = find_cut(keypts, ctr);
imgnew = img(ctr(1)-cut:ctr(1)+cut, ctr(2)-cut:ctr(2)+cut);
% figure;
% imshow(imgnew);title('Image to calculate tilt angle');

%%%%%%%%%%%% determine whether the picture is valid %%%%%%%%%%%%%%%%%%%%%
if ceil(cut)>= win_size
    valid = true;
%     disp('cut length for tilt angle calculation');
%     disp(cut*2+1);
else
    valid = false;
    disp('This is an invalid image!')
end


end
