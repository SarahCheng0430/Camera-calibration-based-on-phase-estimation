function img_pattern = generate_pattern(M, N)
% Usage: generate customized pattern to calibrate
    if  M/2 == round(M/2) || N/2 == round(N/2)
        disp('M and N should be both odd numbers.')
        img_pattern = 0;
        return;
    end
    
    img_pattern = zeros(M * 40, N * 40);
    for ip = 1:1:M
        for jp = 1:1:N
            img_pattern(40*(ip-1)+7:40*(ip-1)+34, 40*(jp-1)+7:40*(jp-1)+34) = 255;
        end
    end
    
    ic = round((M+1)/2);
    jc = round((N+1)/2);
    img_pattern(40*(ic-1)+17:40*(ic-1)+24, 40*(jc-1)+17:40*(jc-1)+24) = 0;
    figure;
    imshow(img_pattern);
end

