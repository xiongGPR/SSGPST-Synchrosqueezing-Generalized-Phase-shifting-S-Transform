function [st,fs,dst] = S_transfom(timeseries,dt)
%  ST
% output:
%   st  ST
%   fs  frequency
%   dst  Differential of ST
% input:
%   timeseries  signal
%   dt  sampling interval
%------------------------------------------------------------------------
%    Authors: Xiong Hongqiang
%    2023/4/1
%---------------------------------------------------------------------------------
timeseries = timeseries(:); % Turn into column vector;
minfreq = 1;     
maxfreq = fix(length(timeseries)/2);
samplingrate=1;
freqsamplingrate=1;
factor = 1;              
%% 
t = (0:length(timeseries)-1)*samplingrate;                      
spe_nelements =ceil((maxfreq - minfreq+1)/freqsamplingrate)   ; 
f = (minfreq + [0:spe_nelements-1]*freqsamplingrate)/(samplingrate*length(timeseries)); 
%% 
n=length(timeseries);
%% Calculate ST and Differential of ST
xMean = mean(timeseries);
timeseries = timeseries-xMean;       
xi = zeros(1, n);                       
xi(1:n/2+1) = 2*pi/n*[0:n/2];           
xi(n/2+2:end) = 2*pi/n*[-n/2+1:-1];     
% the core
vector_fft=fft(timeseries);        %fft
vector_fft=[vector_fft,vector_fft]; %pad 2D
st=zeros(length(f),n);  
if minfreq == 0
    st(1,:) = mean(timeseries)*(1&[1:1:n]);  %f = 0
else
    psih=vector_fft(minfreq+1:minfreq+n).*g_window(n,minfreq,factor);
    st(1,:)=ifft(psih);
    dpsih=(i*xi/dt) .* psih;
    dst(1,:)=ifft(dpsih);
end
%the actual calculation of the ST
% Start loop to increment the frequency point
for banana=freqsamplingrate:freqsamplingrate:(maxfreq-minfreq)  
    psih=vector_fft(minfreq+banana+1:minfreq+banana+n).*g_window(n,minfreq+banana,factor);
    st(banana/freqsamplingrate+1,:)=ifft(psih);
    dpsih=(i*xi/dt) .* psih;
    dst(banana/freqsamplingrate+1,:)=ifft(dpsih);
end   
st1 = abs(st);
st2 = angle(st);
fs = f/dt;

function gauss=g_window(length,freq,factor)
%% new—by xiong
endvector = zeros(1,length);
vector(1:floor(length/2)+1)=[0:floor(length/2)];
vector(floor(length/2)+1:length)=[length-floor(length/2):-1:1];
vector=vector.^2;
vector=vector*(-factor*2*pi^2/freq^2);
gauss=exp(vector);
end
end     
