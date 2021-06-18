%--------------------------------------------------------------------------
%* @author     孟子喻
%* @time       2021.4.20
%* @dependence Communication Toolbox
%*             Image Processing Toolbox
%* @file       audio_FM.m
%*             胡桃.wav
%--------------------------------------------------------------------------


[y, Fs] = audioread("胡桃.wav");
% sound(y,Fs);
y = y';
dt = 1/Fs;
time_during = length(y) * dt; 
t = 0:dt:time_during-dt; 
prinmary_signal = y;
primary_max_A = max(max(prinmary_signal), -min(prinmary_signal));
kf=0.2;             
fc = 200000;        


%---------------------------信号调制----------------------------------------
% 先把每个dt对应的积分求出来备用
integradl_siginal = [];
% 当矩阵元素过多时，最好预分配内存，比如全赋为0等，否则速度会慢很多
start = 1;
step  = 1;
for iter=start:step:length(y)
    integradl_siginal = [integradl_siginal, prinmary_signal(1,iter)];
    % 每个时间内integradl_siginal对应当前时间的积分  
end
modulated_signal = cos(2*pi*fc*t + kf*integradl_siginal);

%--------------------------信号解调-----------------------------------------
% 解调部分参考了网上的代码
% https://blog.csdn.net/mddh_123/article/details/108700454
modulated_signal = diff(modulated_signal);
demodulated_signal = [0 abs(hilbert(modulated_signal))];
demodulated_signal = demodulated_signal - mean(demodulated_signal);% 滤除直流分量


modulated_signal = cos(2*pi*fc*t + kf*integradl_siginal);
% diff处理后modulated_signal会变短一位，所以这里重新算一下，把plot都放在一起
subplot(411);
plot(t,prinmary_signal);
title('基带信号');
xlabel('时间');
ylabel('幅度');

subplot(412);
plot(t,cos(2*pi*fc*t));
title('载波信号');
xlabel('时间');
ylabel('幅度');

subplot(413);
plot(t,modulated_signal);
title('FM调制信号');
xlabel('时间');
ylabel('幅度');
 
subplot(414)
plot(t,demodulated_signal*20);
ylim([-0.5, 0.5]);
title('FM解调信号');
xlabel('时间');
ylabel('幅度');

sound(demodulated_signal,Fs)

