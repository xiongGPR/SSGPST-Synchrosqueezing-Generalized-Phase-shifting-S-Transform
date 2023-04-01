function [Tx, fs,Wx] = synsq_st_fw(x, dt)
%  SS-ST
% output:
%   Tx  SSST
%   fs  frequency
%   Wx  ST
% input:
%   x  signal
%   dt  sampling interval
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------
%% 参数检查
if nargin<2, error('Too few input arguments'); end
%% Calculate ST
[Wx,fs,dWx] = S_transfom(x,dt);  
w= phase_st(Wx,fs,dWx);%Calculate instantaneous frequency
%% Calculate the synchrosqueezed frequency decomposition 
	[na, N] = size(Wx);
    f = 1/size(Wx,2):1/size(Wx,2):size(Wx,1)/size(Wx,2);
    scaleterm = f;
	%incorporate threshold by zeroing out Inf values, so they get ignored below
	Wx(isinf(w)) = 0;      
	Tx = zeros(length(fs),size(Wx,2));
    for b=1:N
        for ai=1:length(fs)
            % Find w_l nearest to w(a_i,b)
            [V,k] = min(abs(w(ai,b)-fs));   
             Tx(k, b) = Tx(k, b) + abs(Wx(ai, b)) /scaleterm(ai);
        end
    end
end


    




