function [rx] = synsq_Fgpst_iw(Tx,A)
% Fast inverse transformation of SS-GPST
% output:
% rx  Reconstructed signal
%
% input:
%    Tx  SS-GPST 
%    A   Generalized parameter
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------
if nargin<2, A = 1/sqrt(3); end
f = 1/size(Tx,2);
rx = 2*sqrt(2*pi)/A*f*sum(real(Tx),1).';
end

 




