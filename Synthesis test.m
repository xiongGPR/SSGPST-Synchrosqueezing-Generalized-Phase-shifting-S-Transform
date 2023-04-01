%% Synthetic signals in the paper
N =1024;%1024个采样点
dt = 10/(N-1);%采样时间
t = 0:dt:10;
h1 = (2+0.2.*cos(t)).*cos(2.*pi.*(3.*t+0.6.*cos(t)));
h2 = 0.7.*(1+0.3.*cos(2.*t)).*exp(-t/20).*cos(2.*pi.*(5.*t+0.6.*t.^2+0.3.*sin(t)));
h3 = 0.4.*cos(2*pi*9.*t);
signal = h1+h2+h3; %信号
% signal = awgn(signal,20,'measured');%加入高斯白噪声
figure(1);
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nexttile;plot(t,signal,'k');ylabel('Amplitude');set(gca,'YLim',[-5 5],'FontName','Times New Roman','FontSize',10); 
nexttile;plot(t,h1,'k');ylabel('Amplitude');set(gca,'YLim',[-3 3],'FontName','Times New Roman','FontSize',10);
nexttile;plot(t,h2,'k');ylabel('Amplitude');set(gca,'YLim',[-2 2],'FontName','Times New Roman','FontSize',10);
nexttile;plot(t,h3,'k');ylabel('Amplitude');set(gca,'YLim',[-2 2],'FontName','Times New Roman','FontSize',10);xlabel('Time (s)');
set(gca,'color','white'); 
%% Time-frequency spectrum
A1 = 1;%一般为1/sqrt(3)
A2 = 0.68;%一般为1/sqrt(3)（约为0.58），0.68合适
x = signal;
[st,fs,dst] = S_transfom(x,dt);  %ST结果
[pgst1,fs,dpgst1,f] = GPST_fw(x,dt,A1);  %PGST结果
[pgst2,fs,dpgst2,f] = GPST_fw(x,dt,A2);  %PGST结果
st_out = abs(st);
pgst1_out = abs(pgst1);
pgst2_out = abs(pgst2);
color = jet;
 %% SS -Time-frequency---spectrum
x = signal;
[ssst, fs,st] = synsq_st_fw(x, dt);%SS-ST
[sspgst1, fs,pgst1] = synsq_gpst_fw(x, dt,A1);%SS-GPST,A1
[sspgst2, fs,pgst2] = synsq_gpst_fw(x, dt,A2);%SS-GPST,A2
 ssst_out = abs(ssst);
 sspgst1_out = abs(sspgst1);
 sspgst2_out = abs(sspgst2);
   %phase spectrum
%  ssst_out = angle(ssst);
%  sspgst1_out = angle(sspgst1);
%  sspgst2_out = angle(sspgst2);
 % 归一化
%注意：只有按照熊的定义才能精确重构，但是在做Fig5时，需要改成黄的那种样子，看起来图的对比更强
 %% 算法测试，成图部分
 f = 20;%指定20Hz内的信号,重要输入参数1
 [~,k] = min(abs(fs-f));%计算频率点长度
 figure(3);
 tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
 %SSST 成图
 nexttile;imagesc( t,fs(1:k+1),ssst_out(1:k+1,:));%绘图语句
ylabel('Frequency (Hz)');xlabel({'Time (s)'});%轴标题
%  annotation('textbox',[.82 .82 .1 .1],'String','SSST','FitBoxToText','on');%注释部分
 map = color; map = map(end:-1:1,:); colormap(map); colorbar;axis xy;%灰度，并调整坐标轴
 set(gca,'FontName','Times New Roman','FontSize',10);
 %SSGST 成图
nexttile;imagesc( t,fs(1:k+1),sspgst1_out(1:k+1,:));%绘图语句
ylabel('Frequency (Hz)');xlabel({'Time (s)'});%轴标题
%  annotation('textbox',[.374 .348 .1 .1],'String','SSGST','FitBoxToText','on');%注释部分
 map = color; map = map(end:-1:1,:); colormap(map); colorbar;axis xy;%灰度，并调整坐标轴
 set(gca,'FontName','Times New Roman','FontSize',10);
 %SSPGST 成图
nexttile;imagesc( t,fs(1:k+1),sspgst2_out(1:k+1,:));%绘图语句
ylabel('Frequency (Hz)');xlabel({'Time (s)'});%轴标题
%  annotation('textbox',[.805 .348 .1 .1],'String','SSPGST','FitBoxToText','on');%注释部分
 map = color; map = map(end:-1:1,:); colormap(map); colorbar;axis xy;%灰度，并调整坐标轴
 set(gca,'FontName','Times New Roman','FontSize',10);
 
 %% 重构误差分析


rx_ssst = synsq_st_iw(ssst, st);%ssst重构，基本正确
rx_ssst_out = abs(rx_ssst(:,1)'-signal);
ssst_mean = mean(rx_ssst_out);%均值
ssst_max = max(rx_ssst_out);%最大值
ssst_mse = immse(rx_ssst(:,1)',signal);%均方误差

rx_sspgst1 = synsq_gpst_iw(sspgst1,A1);%sspgst1重构，正确
rx_sspgst1_out = abs(rx_sspgst1(:,1)'-signal);
sspgst1_mean = mean(rx_sspgst1_out);%均值
sspgst1_max = max(rx_sspgst1_out);%最大值
sspgst1_mse = immse(rx_sspgst1(:,1)',signal);%均方误差

rx_sspgst2 = synsq_gpst_iw(sspgst2,A2);%sspgst2重构，正确
rx_sspgst2_out = abs(rx_sspgst2(:,1)'-signal);
sspgst2_mean = mean(rx_sspgst2_out);%均值
sspgst2_max = max(rx_sspgst2_out);%最大值
sspgst2_mse = immse(rx_sspgst2(:,1)',signal);%均方误差

rx_Fsspgst1 = synsq_Fgpst_iw(sspgst1,A1);%Fsspgst1重构，正确
rx_Fsspgst1_out = abs(rx_Fsspgst1'-signal);
Fsspgst1_mean = mean(rx_Fsspgst1_out);%均值
Fsspgst1_max = max(rx_Fsspgst1_out);%最大值
Fsspgst1_mse = immse(rx_Fsspgst1',signal);%均方误差

rx_Fsspgst2 = synsq_Fgpst_iw(sspgst2,A2);%Fsspgst2重构，正确
rx_Fsspgst2_out = abs(rx_Fsspgst2'-signal);
Fsspgst2_mean = mean(rx_Fsspgst2_out);%均值
Fsspgst2_max = max(rx_Fsspgst2_out);%最大值
Fsspgst2_mse = immse(rx_Fsspgst2',signal);%均方误差



figure(4);%重构信号图
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
plot(t,signal,'--',t,rx_ssst(:,1),'g',t,rx_sspgst1(:,1),'b',t,rx_sspgst2(:,1),'r');
legend('s(t)','SS-ST','SS-GPST(A=1)','SS-GPST(A=0.68)','NumColumns',3,'FontSize',8);ylabel('Amplitude');xlabel({'Time (s)'});%轴标题
set(gca,'YLim',[-8 8],'FontName','Times New Roman','FontSize',10); 
nexttile
plot(t,rx_ssst_out,'g',t,rx_sspgst1_out,'b',t,rx_sspgst2_out,'r');
legend('SST','SS-ST','SS-GPST(A=1)','SS-GPST(A=0.68)','NumColumns',2,'FontSize',8);ylabel('Amplitude');xlabel({'Time (s)'});%轴标题 
set(gca,'YLim',[-1 1],'FontName','Times New Roman','FontSize',10); 



