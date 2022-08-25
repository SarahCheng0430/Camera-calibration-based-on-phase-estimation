%% -------------------- PE2D_length3.m -----------------------
% Usage:Calculate PE2D for sub-image, with a window sliding along long
%       edge, to determine the world coorinate      
%

function [phase1,phase2]=PE2D_length3(N_in,data,col_num,row_start,row_end,col_start,col_end)
    i=0;
    for  k=col_start:N_in-col_end+1 
        data_im=data(row_start:row_end,k:k+col_num-1);
        i=i+1;
        [N,M]=size(data_im);
        midpoint=ceil([(N+1)/2 (M+1)/2]);
        win2=@hann;
        win1=@hann;
        win=window2(N,M,win1,win2);
        L_square = 10;
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%calculate the first sub-image
        if(i==1)
            [phase1,phase2,u1,v1,u2,v2] = PE2D(data_im);
            phase1(i)=phase1;
            phase2(i)=phase2;
            u1(i) = u1;v1(i) = v1;u2(i) = u2; v2(i) = v2; 
          
            row1(i) = u1(i)+midpoint(1);
            col1(i) = v1(i)+midpoint(2);
            row2(i) = u2(i)+midpoint(1);
            col2(i) = v2(i)+midpoint(2);
        end

        %%calculate the next sub-images
        if(i~=1)
            Spectrum=fftshift(fft2(data_im.*win));
            S_value=abs(Spectrum);
            
        %search the maximum points in spectrum
            SearchSquare_top_1(i) = max((row1(i-1)-L_square),1);
            SearchSquare_bottom_1(i) = min((row1(i-1)+L_square),N);
            SearchSquare_left_1(i) = max((col1(i-1)-L_square),1);
            SearchSquare_right_1(i) = min((col1(i-1)+L_square),M);

            SearchSquare_top_2(i) = max((row2(i-1)-L_square),1);
            SearchSquare_bottom_2(i) = min((row2(i-1)+L_square),N);
            SearchSquare_left_2(i) = max((col2(i-1)-L_square),1);
            SearchSquare_right_2(i) = min((col2(i-1)+L_square),M);

            search_mat_1_D1=S_value(SearchSquare_top_1(i):SearchSquare_bottom_1(i),SearchSquare_left_1(i):SearchSquare_right_1(i));
            search_mat_1_D2=S_value(SearchSquare_top_2(i):SearchSquare_bottom_2(i),SearchSquare_left_2(i):SearchSquare_right_2(i));
            
          %---------------------------------------------------------------
            [row1_s(i), col1_s(i)] = find(search_mat_1_D1==max(max(search_mat_1_D1)));
            row1(i)=row1_s(i)+SearchSquare_top_1(i)-1;
            col1(i)=col1_s(i)+SearchSquare_left_1(i)-1;
          
            u1(i)=row1(i)-midpoint(1);
            v1(i)=col1(i)-midpoint(2);
            point1=2;
            k1(i)=(u1(i)/N)/(v1(i)/M);
            theta1(i)=atan(1/k1(i))/pi*180;

          %----------------------------------------------------------------
            [row2_s(i), col2_s(i)] = find(search_mat_1_D2==max(max(search_mat_1_D2)));
            row2(i)=row2_s(i)+SearchSquare_top_2(i)-1;
            col2(i)=col2_s(i)+SearchSquare_left_2(i)-1;
        
            u2(i)=row2(i)-midpoint(1);
            v2(i)=col2(i)-midpoint(2);
            point2=2;
          if u2(i)>0
            u2(i)=row2(i)-midpoint(1);
            v2(i)=col2(i)-midpoint(2);
            point2=1;
          end
            k2(i)=(u2(i)/N)/(v2(i)/M);
            theta2(i)=atan(1/k2(i))/pi*180;
    
           % obtain the phases of next sub-images
            phase1(i)=angle(Spectrum(row1(i),col1(i)))+pi*u1(i)*(N-1)/N+pi*v1(i)*(M-1)/M;
            phase1(i)=mod(phase1(i),2*pi);
            phase2(i)=angle(Spectrum(row2(i),col2(i)))+pi*u2(i)*(N-1)/N+pi*v2(i)*(M-1)/M;
            phase2(i)=mod(phase2(i),2*pi);
           
        end

    end % k end
    
end % function end
