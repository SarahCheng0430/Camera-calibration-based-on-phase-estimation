function err_r = main_re_app(Data)
%  - Usage: calculate mean reprojection error
%  - Plot the figures of MRE of each input image and RMS of control points
clc;
P = Data.p; 
Q = Data.q;
In = Data.In; 
k = Data.k; 
Rm = Data.Rm;
option = Data.opt;

du = In(1,3);
dv = In(2,3);
[~,~,num] = size(Rm);

figure;
hold on;
err_t1 = [];
err_t2 = [];
for i = 1:num    
    [ori_pts, XY] = calc_projection(P(:,:,i), In, Rm(:,:,i));
%     option1
    ori_pts = undistort1(ori_pts, k, XY, du, dv);    
    tar_pts = Q(:,:,i)';
%     option2
%   tar_pts = undistort2(tar_pts, k, T, du, dv);
    [err1, err2, x, y] = calc_re(ori_pts, tar_pts);
    err_t1 = [err_t1, err1];
    err_t2 = [err_t2, err2];
    scatter(x,y,'filled');    
end
title(option + " - number of figures: " + string(num));

figure;
b=bar(err_t1);
grid on;
legend('reprojection error');
xlabel('picture');
ylabel('reprojection error');
title(option + " - number of figures: " + string(num));
disp(option + " - mean reprojection error: " + string(sum(err_t1) / num)); %mean reprojection error
disp(option + " - root mean square error: " + string(sqrt(sum(err_t2) / num)));%root mean square error
err_r = sum(err_t1) / num;

    function [projected_pts, XY] = calc_projection(pts, In, Ex)
        pts(3,:) = 1;
        XY = Ex * pts;
        projected_pts = In * XY;
        projected_pts = projected_pts(1:2,:) ./ projected_pts(3,:);
        projected_pts = projected_pts';
        XY = XY(1:2,:) ./ XY(3,:);
        XY = XY';
    end

    function [err1, err2 , x, y] = calc_re(ori_pts, tar_pts) 
        err_pts = ori_pts - tar_pts;
        err = err_pts(:,1) .* err_pts(:,1) + err_pts(:,2) .* err_pts(:,2);
        err1 = sum(sqrt(err))/length(err);
        err2 = sum(err)/length(err);
        x = err_pts(:,1);
        y = err_pts(:,2);
    end

    function corrected = undistort1(pts, k, XY, du, dv)
        k1 = k(1); k2 = k(2);
        r2 = XY(:,1) .* XY(:,1) + XY(:,2) .* XY(:,2);
        X = pts(:,1) + (pts(:,1) - du) .* (k1 * r2 + k2 * r2 .* r2);
        Y = pts(:,2) + (pts(:,2) - dv) .* (k1 * r2 + k2 * r2 .* r2);
        corrected = [X,Y];
    end

end

