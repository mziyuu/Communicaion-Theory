%--------------------------------------------------------------------------
%* @author     孟子喻
%* @time       2021.4.20
%* @dependence Communication Toolbox
%*             Image Processing Toolbox
%* @file       audio_DSB.m
%*             胡桃.wav
%* @reference  lowpass用法参考了https://blog.csdn.net/weixin_43870101/article/details/106888716
%*             在从网络上搜索到一些方法以后，我觉得很多方法存在的问题是没能好好
%*             利用已有的函数，自己造轮子比较多，所以这篇代码几乎没有参考其他
%*             资料，单纯使用课本公式搭配官方的各种函数实现。
%--------------------------------------------------------------------------




%---------------------读取wav音频文件---------------------------------------
[y, Fs] = audioread("胡桃.wav");
% sound(y,Fs);
y = y';
dt = 1/Fs;
time_during = length(y) * dt; 
t = 0:dt:time_during-dt; 
% 这里按理说不需要减去dt，但是可能是因为精确度的问题
% t多了一个维度（应是443392.实为443393），所以要减去一个保证维度一致



%--------------------载波调制-----------------------------------------------
fc = 200000;  %不能使用Fs整数倍，一般为固定值
prinmary_signal = y;
modulated_signal =  prinmary_signal.*cos(2*pi*fc*t);



%--------------------相干滤波-----------------------------------------------
demodulated_signal_before_lpf = (modulated_signal.*cos(2*pi*fc*t)-0.5)*2;
demodulated_signal_after_lpf = lowpass(demodulated_signal_before_lpf, Fs/2, fc);

%在信号x中加入高斯白噪声，信噪比SNR以dB为单位
modulated_signal_before_lpf_noise = awgn(modulated_signal, 20);
demodulated_signal_before_lpf_noise = (modulated_signal_before_lpf_noise.*cos(2*pi*fc*t)-0.5)*2;
demodulated_signal_after_lpf_noise = lowpass(demodulated_signal_before_lpf, Fs/2, fc);



%----------------绘图区-----------------------------------------------------


figure
subplot(711)
plot(t, prinmary_signal)
xlabel("时间")
ylabel("幅度")
title("基带信号")

subplot(712)
plot(t, cos(2*pi*fc*t))
xlabel("时间")
ylabel("幅度")
title("载波信号")

subplot(713)
plot(t, modulated_signal)
xlabel("时间")
ylabel("幅度")
title("调制信号")

subplot(714)
plot(t, demodulated_signal_before_lpf)
xlabel("时间")
ylabel("幅度")
title("解调信号-滤波前")

subplot(715)
plot(t, demodulated_signal_after_lpf)
xlabel("时间")
ylabel("幅度")
title("解调信号-滤波后")

subplot(716)
plot(t, demodulated_signal_before_lpf_noise)
xlabel("时间")
ylabel("幅度")
title("解调信号-滤波前-加高斯噪声")

subplot(717)
plot(t, demodulated_signal_after_lpf_noise)
xlabel("时间")
ylabel("幅度")
title("解调信号-滤波后-加高斯噪声")

%---听声音对比----
sound(demodulated_signal_after_lpf, Fs)
%sound(demodulated_signal_after_lpf_noise, Fs)

y = modulated_signal;
z = abs(fft(y));
m = z(1:length(z)/2);
n = length(y);          
f = (1:n/2)*(Fs/n); 
figure
plot(f,m)
xlabel('Frequency')
ylabel('Power')
%参考https://blog.csdn.net/Yujie_Wang/article/details/112621679