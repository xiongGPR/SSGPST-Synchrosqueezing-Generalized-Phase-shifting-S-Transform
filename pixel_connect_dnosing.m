
function dnosing = pixel_connect_dnosing(f_sst)
%  Pixel connectivity filtering
% output:
%   dnosing  filter matrix 
% input:
%   f_sst  Time-frequency spectrum
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------



     %step 1    Hard threshold filtering
    [na, n] = size(f_sst);
    Wx_fine = abs(f_sst);
    lamda = sqrt(2*log(n)) * mad(Wx_fine) * 1.4286;
    f_sst(abs(f_sst)<lamda)=0;   
%     f_sst(abs(f_sst)~=0)=1;
%     dnosing = f_sst;
    
    
    %step 1      Connectivity filtering       
    f_sst_con = abs(f_sst);
    [f_sst_con,L] = bwlabel(f_sst_con);

            m=zeros(L,1);
            for i= 1:L                    
                m(i) = sum(f_sst_con(:)==i);  
            end
            % Set to 0 if less than the threshold
            mm = mad(m);
            mmm = sqrt(2*log(L)) * mm * 0.2;   %  the threshold
            % if 
            for i= 1:L
                if m(i)<mmm
                    f_sst_con(f_sst_con == i) = 0;
                end
            end
                f_sst_con(f_sst_con>0) = 1; 
            dnosing =  f_sst_con;  %Return filter matrix 

end


