function [rx] = synsq_gpst_iw(Tx,A,Cs,freqband)
% Fast inverse transformation of SS-GPST
% output:
% rx  Reconstructed signal and error
%
%
% input:
%    Tx  SS-GPST 
%    A   Generalized parameter
%    Cs (optional) curve centerpoints  
%    freqband (optional) curve bands       
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------
if nargin<4, freqband = size(Tx,1); end
if nargin<3, Cs = ones(size(Tx,2),1); end 
if nargin<2, A = 1/sqrt(3); end
  C1 = quadgk(@(w) 0.5*exp(-0.5/A^2*(w-2*pi).^2)./w, 0, Inf); %tolerance factor
%% 
f = 1/size(Tx,2):1/size(Tx,2):size(Tx,1)/size(Tx,2);
rx = zeros(size(Cs,1),size(Cs,2)+1);  
TxMask = zeros(size(Tx));
TxRemainder = Tx;    
	for n=[1:size(Cs,2)]   
		TxMask = zeros(size(Tx));
		UpperCs=min(max(Cs(:,n)+freqband(:,n),1),length(f));   
		LowerCs=min(max(Cs(:,n)-freqband(:,n),1),length(f));  
		%Cs==0 corresponds to no curve at that time, so this removes such points from the inversion
		UpperCs(find(Cs(:,n)<1))=1;
		LowerCs(find(Cs(:,n)<1))=2;
		for m=[1:size(Tx,2)]  
			TxMask(LowerCs(m):UpperCs(m),m) = Tx(LowerCs(m):UpperCs(m), m);
			TxRemainder(LowerCs(m):UpperCs(m),m) = 0;            
        end
		rx(:,n) = 1*(f(2)-f(1))/C1*sum(real(TxMask),1).';  
	end
	rx(:,n+1) = 1*(f(2)-f(1))/C1*sum(real(TxRemainder),1).'; 
end

 




