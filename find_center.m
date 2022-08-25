function [u,v,imgcut] = find_center(img)
% - Usage: determination of the world coordinate origin in the image plane
% - Input: img (remapped image)
% - Output: [u(h),v(w)] of the world coordinate origin Ow
    
%% connected region detection for each grid
    function [keypts, bd] = find_bd(Is)
        I = edge(Is, 'log');
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
        [h,w] = size(imgnew);
        m1 = h+w;
        m2 = 0;
        n1 = h-1;
        n2 = 1-w;
        keypts = zeros(4,2);
        for x1 = 1:h
            for y1 = 1:w
                if imgnew(x1,y1) == true
                    sum_m = x1 + y1;
                    sum_n = x1 - y1;
                    if sum_m < m1
                        keypts(1,:) = [x1,y1];
                        m1 = sum_m;
                    end
                    if sum_m > m2
                        keypts(2,:) = [x1,y1];
                        m2 = sum_m;
                    end
                    if sum_n < n1
                        keypts(3,:) = [x1,y1];
                        n1 = sum_n;
                    end
                    if sum_n > n2
                        keypts(4,:) = [x1,y1];
                        n2 = sum_n;
                    end
                end
            end
        end
        left = min([keypts(1,2),keypts(4,2)]);
        right = max([keypts(2,2),keypts(3,2)]);
        up = min([keypts(1,1),keypts(3,1)]);
        down = max([keypts(2,1),keypts(4,1)]);
        bd = round([up, down, left, right]);       
    end
    
%% find the pixel coordinate of the pattern center
    function [II, ctr] = find_center(Ic, ii, jj)
        lst = [];
        ls = [ii,jj];
        Ic(ii,jj) = false;
        
        while true
            [N,~] = size(ls);
            for n = 1:N
                lst = [lst; ls(n,:)];
                iii = ls(n,1);
                jjj = ls(n,2);       
                if Ic(iii-1,jjj-1) == true
                    ls = [ls; iii-1,jjj-1];
                    Ic(iii-1,jjj-1) = false;
                end
                if Ic(iii-1,jjj) == true
                    ls = [ls; iii-1,jjj];
                    Ic(iii-1,jjj) = false;
                end
                if Ic(iii-1,jjj+1) == true
                    ls = [ls; iii-1,jjj+1];
                    Ic(iii-1,jjj+1) = false;
                end
                if Ic(iii,jjj-1) == true
                    ls = [ls; iii,jjj-1];
                    Ic(iii,jjj-1) = false;
                end
                if Ic(iii,jjj+1) == true
                    ls = [ls; iii,jjj+1];
                    Ic(iii,jjj+1) = false;
                end
                if Ic(iii+1,jjj-1) == true
                    ls = [ls; iii+1,jjj-1];
                    Ic(iii+1,jjj-1) = false;
                end
                if Ic(iii+1,jjj) == true
                    ls = [ls; iii+1,jjj];
                    Ic(iii+1,jjj) = false;
                end
                if Ic(iii+1,jjj+1) == true
                    ls = [ls; iii+1,jjj+1];
                    Ic(iii+1,jjj+1) = false;
                end               
            end
            [N1,~] = size(ls);
            if N1 == N
                break;
            else
                ls = ls(N+1:N1,:);
            end
        end
        [Np,~] = size(lst);
        if abs(max(lst(:,1))-min(lst(:,1)) - (max(lst(:,2))-min(lst(:,2)))) <= 4 && Np >= 15
%            ctr = [mean(lst(:,1)), mean(lst(:,2))];
             ctr = [mean([max(lst(:,1)), min(lst(:,1))]), mean([max(lst(:,2)), min(lst(:,2))])];
        else
            ctr = [0,0];
        end
        II = Ic;
    end

    function cpt = closest_pt(pts, pt)
       [NE,~] = size(pt);
       if NE > 1
           pt = pt(1,:);
       end
        err = pts - pt;
        norm = err(:,1).* err(:,1) + err(:,2).* err(:,2);
        cpt = pts(norm == min(norm),:);    
        
        %if there are more than one suited points, choose the first point.
        [Nss,~] = size(cpt);
        if Nss > 1
            cpt = cpt(1,:);
        end
    end

    function npt = find_next(pts, pt, inc)
        is = 0;
        while true
            xp = closest_pt(pts, pt + is .* inc);
            if is >= 20
                npt = [0,0];
                break;
            end
            if all(xp == pt) == true
                is = is + 1;
            else
                npt = xp;
                break;
            end
        end
    end

%% crop the remapped image for world coordinate calculation
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
        while cut <= 400
            l = (cut+19)/2;
            j1 = judge(mainpts, center+[l,l]);
            j2 = judge(mainpts, center+[-l,l]);
            j3 = judge(mainpts, center+[l,-l]);
            j4 = judge(mainpts, center+[-l,-l]);
            if j1 && j2 && j3 && j4
                cut = cut + 20;
            else
                break;
            end
        end
        cut = (cut-1)/2;
    end

%% main 
    [keypts, bd] = find_bd(img);
    It = edge(img, 'log');
%     figure;
%     imshow(It);
    
    centers = [];% centers of all grids
%     disp(bd);
    for i = bd(1):bd(2)
        for j = bd(3):bd(4)
            if It(i,j) == true
                [It, ct] = find_center(It, i, j);
                if ct(1) > 0
                    centers = [centers; ct];%ct--count
                end
            end
        end
    end
    
    Jud = centers(:,1) + centers(:,2);
    start = centers(Jud == min(Jud),:);
    H = 1; vec = start(1,:);
    while true
        vec = find_next(centers, vec, [1,0]);
        if vec(1) > 0
            H = H + 1;
        else
            break;
        end
    end
    W = 1; hor = start(1,:);
    while true
        hor = find_next(centers, hor, [0,1]);
        if hor(1) > 0
            W = W + 1;
        else
            break;
        end
    end
    
    for step = 1:round((H-1)/2)
        start = find_next(centers, start, [1,0]);
    end
    for step = 1:round((W-1)/2)
        start = find_next(centers, start, [0,1]);
    end
    %img(round(start(1)),round(start(2))) = 0;
    
    % function output 
    u = start(1);v = start(2);
    cut = find_cut(keypts, start);
    
%     disp('pattern center');
%     disp(u);disp(v);
%     disp('cut length for point calculation');
%     disp(cut*2+1);
    imgcut = img(u-cut:u+cut, v-cut:v+cut);
%     figure;
%     imshow(imgcut);%
    
end

