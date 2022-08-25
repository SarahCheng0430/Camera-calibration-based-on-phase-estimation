function ptss = extract_pt(imgss)
%% connected region detection
    function [bd] = find_bd(Is)
        I = edge(Is, 'log');
        kernel = ones(10,10); 
        I_s = imdilate(I, kernel);
        I_s = imerode(I_s, kernel);
        I_s = imerode(I_s, kernel);
        I_s = imdilate(I_s, kernel);
        imLabel = bwlabel(I_s);
        stats = regionprops(imLabel,'Area');   % find the size of each connected area
        area = cat(1,stats.Area);
        index = find(area == max(area));  % find the index of the largest connected area      
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
        bd = [up, down, left, right];       
    end
    
%% centriod calculation for all grids
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
        if abs(max(lst(:,1))-min(lst(:,1)) - (max(lst(:,2))-min(lst(:,2)))) <= 4 && Np > 10
            ctr = fit_square(lst);
        else
            ctr = [0,0];
        end
        II = Ic;
    end

%% 
    function cpt = closest_pt(pts, pt)
        err = pts - pt;
        norm = err(:,1).* err(:,1) + err(:,2).* err(:,2);
        cpt = pts(norm == min(norm),:);        
        [Nss, ~] = size(cpt);
        if Nss > 1
            cpt = cpt(1,:);
        end
    end

    function npt = find_next(pts, pt, inc)
        if length(size(pt)) == 3
            pt = [pt(1,1,1),pt(1,1,2)];
        end
        is = 0;
        while true
            xp = closest_pt(pts, pt + is .* inc);
            if size(xp,1) == 2
            end
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

%% main
img = imgss;
bd = find_bd(img);
It = edge(img, 'log');
%imshow(It);

% find center points
centers = [];
for i = bd(1):bd(2)
    for j = bd(3):bd(4)
        if It(i,j) == true
            [It, ct] = find_center(It, i, j);
            if ct(1) > 0
                centers = [centers; ct];
            end
        end
    end
end

[ts,~] = size(centers);
for thes = 1:ts
    It(round(centers(thes,1)),round(centers(thes,2))) = true;
end
%imshow(It);

Jud = centers(:,1) + centers(:,2);
start = centers(Jud == min(Jud),:);

% find size
H = 1; vec = start;
while true
    vec = find_next(centers, vec, [1,0]);
    if vec(1) > 0
        H = H + 1;
    else
        break;
    end
end
W = 1; hor = start;
while true
    hor = find_next(centers, hor, [0,1]);
    if hor(1) > 0
        W = W + 1;
    else
        break;
    end
end

% extract pts
pt_lattice = zeros(H,W,2);
pt_lattice(1,1,:) = start;
for ip = 2:W
    pt_lattice(1,ip,:) = find_next(centers, pt_lattice(1,ip-1,:), [0,1]);
end
for ip1 = 2:H
    pt_lattice(ip1,1,:) = find_next(centers, pt_lattice(ip1-1,1,:), [1,0]);
    for ip2 = 2:W
        pt_lattice(ip1,ip2,:) = find_next(centers, pt_lattice(ip1,ip2-1,:), [0,1]);
    end
end

ptss = pt_lattice;

end