function [rx] = synsq_st_iw(Tx,Wx, Cs,freqband)
% Fast inverse transformation of SS-GPST
% output:
% rx  Reconstructed signal and error
%
% input:
%    Tx  SS-GPST 
%    Cs (optional) curve centerpoints  
%    freqband (optional) curve bands       
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------
if nargin<4, freqband = size(Tx,1); end
if nargin<3, Cs = ones(size(Tx,2),1); end 
C1 = quadgk(@(x) exp(-1/2*(x-2*pi).^2)./x, 0, Inf);
C1 = C1*0.5;
C2 = zeros(size(Tx)) ;
t = 0:1:size(Tx,2)-1;
f = 1/size(Tx,2):1/size(Tx,2):size(Tx,1)/size(Tx,2);
% Huang's method
% part1 =1./f.^2;
% part1 = part1'*ones(1,size(Tx,2));
% part2 = angle(Wx)+2*pi*f'*t;
% part2 = exp(i*part2);
% C2 = part1.*part2;
% % Xiong's method
part2 = angle(Wx)+2*pi*f'*t;
C2 = exp(i*part2);

%% 
rx = zeros(size(Cs,1),size(Cs,2)+1);  
TxMask = zeros(size(Tx));
TxRemainder = Tx.*C2;    
	for n=[1:size(Cs,2)]   
		UpperCs=min(max(Cs(:,n)+freqband(:,n),1),length(f));   
		LowerCs=min(max(Cs(:,n)-freqband(:,n),1),length(f));   
		%Cs==0 corresponds to no curve at that time, so this removes such points from the inversion
		UpperCs(find(Cs(:,n)<1))=1;
		LowerCs(find(Cs(:,n)<1))=2;
		for m=[1:size(Tx,2)]  
			TxMask(LowerCs(m):UpperCs(m),m) = Tx(LowerCs(m):UpperCs(m), m).*C2(LowerCs(m):UpperCs(m), m);
			TxRemainder(LowerCs(m):UpperCs(m),m) = 0;
        end    
        rx(:,n) = 1/C1*(f(2)-f(1))*sum(real(TxMask),1).';   
	end
	rx(:,n+1) = 1/C1*(f(2)-f(1))*sum(real(TxRemainder),1).'; 

end



    




