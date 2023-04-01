function w= phase_st(Wx, fs,dWx)
% Calculate instantaneous frequency
%   w       instantaneous frequency
%   Wx      ST
%   fs      frequency
%   dWx     Differentiation of ST
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------

[M,N] = size(Wx);
w = (imag(dWx./Wx/(2*pi)))+fs'*ones(1,N); 
w(abs(Wx)<eps) = Inf;   
end

