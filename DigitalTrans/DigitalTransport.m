%--------------------------------------------------------------------------
% @auther      孟子喻
% @file        transport.m
% @time        2021.5.20
% @dependence  ImageProcessingToolbox   ---awgn()
%              MATLAB >= 2018.a         ---lowpass()
% @reference   MathWorks
%--------------------------------------------------------------------------
clear all
%------使用蒙特卡洛法分析--------------------------------------------------
% 查看多次结果
ask_acc_list = [];
psk_acc_list = [];
snr_list = [];
for snr = 1:1:30
    ask_acc_list_temp = [];
    psk_acc_list_temp = [];  % 临时存储单个准确率，计算出平均值后再放入该信噪比下的准确率
    for iter=1:10
        [acc_ask,acc_psk] = DigitalTransport(iter,snr);  % MATLAB中多个返回值要用[]
        ask_acc_list_temp = [ask_acc_list_temp, acc_ask];
        psk_acc_list_temp = [psk_acc_list_temp, acc_psk];
    end
    ask_acc_list = [ask_acc_list,mean(ask_acc_list_temp)];
    psk_acc_list = [psk_acc_list,mean(psk_acc_list_temp)];
    
    snr_list = [snr_list, snr];
end
%------绘制不同信噪比下各传输方式准确率------------------------------------
% 2ASK存在一个问题就是无法找到在所有信噪比下都最优的阈值，所以变化明显

figure(4)
plot(snr_list, ask_acc_list, '-')
hold on;
plot(snr_list, psk_acc_list, '-.')
xlabel("信噪比")
ylabel("准确率")
ylim([0,110])
legend("2ASK","2PSK")
title("不同信噪比下各传输方式准确率")

function [acc_ask,acc_psk] = DigitalTransport(test_iter, snr)
% MATLAB中多个返回值要用[]，函数声明和调用时同
str = 'aaabbiieubaddmengziyu';
fprintf("\nTest iter %d \t snr %d\n", test_iter,snr)
%-------Huffman编码--------------------------------------------------------
str_len = length(str);
char_type = unique(str);              %unique函数用于获取字符串中出现过的字符
char_type_num = length(char_type);    %测量有多少种字符

char_type_cell = cell(1, char_type_num);
p = zeros(1, char_type_num);

% 将字符类型存入元胞数组，方便huffmandict调用
for i = 1:char_type_num
    char_type_cell{1,i} = char_type(i);
end

% 计算各字符出现概率
for i = 1:char_type_num
    p(i) = numel(find(str==char_type(i))) / str_len;
end


%---Huffman字典---
dict = huffmandict(char_type_cell, p);    % 返回的dict储存huffman字典（即对应关系）
if test_iter == 1
    for i = 1:length(char_type)
        fprintf('%s: %s', char_type(i),string(dict{i,2}));
            % dict是一个元胞数组cell，第一个元素是字符，第二个元素是对应的Huuffman码
        fprintf("\n")
    end
end

fprintf("发送方字符串\t\t\t\t\t:")
fprintf("%s", str)
fprintf("\n")

%---Huffman码---
huff_code = huffmanenco(str, dict);            % 根据dict和str进行Huffman编码
if test_iter == 1
    fprintf('发送方Huffman编码\t\t\t:');
    for i = 1:str_len
        fprintf('%d', dict{char_type==str(i),2});
    end
    fprintf('\n');
end

%-------Hamming编码--------------------------------------------------------
n = 7;
k = 4;
ham_code = encode(huff_code,n,k,'hamming/binary');
if test_iter == 1
    fprintf("发送方Hamming编码\t\t\t:")
    fprintf("%s", string(ham_code))
    fprintf("\n")
    fprintf("\n")
end
%-------Hamming直接解码-----------------------------------------------------
huff_code = decode(ham_code, n, k, 'hamming/binary');
if test_iter == 1
    fprintf("接收方Hamming译码(参考)\t\t:")
    fprintf("%s", string(huff_code))
    fprintf("\n")
end

%-------Huffman直接解码-----------------------------------------------------
str_decoded = huffmandeco(huff_code, dict); % 根据huff_code和dict还原字符串
if test_iter == 1
    fprintf("接收方Huffman译码(参考)\t\t:")
    fprintf("%s", string(str_decoded))
    fprintf("\n")
end

%------原信号--------------------------------------------------------------
T = 1;          % 一个脉冲周期
sample = 100;   % 一个脉冲周期内的采样点
mt=[];          % 原信号
t=0+(T/sample):(T/sample):length(ham_code);
for i=1:length(ham_code)          
    if(ham_code(i)==1)
        % 一个周期内设置为1
        for j = 1:sample
            mt = [mt, 1];
        end
    else    %ham_code(i)==0
        % 一个周期内设置为0
        for j = 1:sample
            mt = [mt, 0];
        end
    end
end
if test_iter == 1
    figure(1)
    subplot(4,1,1)
    plot(t,mt);
    ylim([-0.5,1.5]);
    title('原信号')
end

%------2ASK编码------------------------------------------------------------
sin_wave = sin(2*pi/T .*t);
mt_ask=[];
for i=1:length(ham_code)
    if(ham_code(i)==1)
        for j = 1:sample
            mt_ask = [mt_ask,sin_wave((i-1)*sample + j)]; %每个区间的开始，所以要-1
        end
    else
        for j = 1:sample
            mt_ask = [mt_ask,0];
        end
    end
end
mt_ask_noise = awgn(mt_ask, snr);  % 加噪声
if test_iter == 1
    subplot(4,1,2)
    plot(t,mt_ask);
    ylim([-1.5,1.5]);
    title('2ASK')
end

%------2ASK译码------------------------------------------------------------
%------使用包络检波译码------
mt_ask_rectify = abs(mt_ask_noise);  %整流
mt_ask_filted  = lowpass(mt_ask_rectify,1/T/100,1/T);

mt_ask_deco = [];
for t_judge = 0.5:T:length(ham_code)-0.5 % t_sampling是抽样判决时刻
    if mt_ask_filted(t_judge*sample)> 0.1
        mt_ask_deco = [mt_ask_deco,1];
    else
        mt_ask_deco = [mt_ask_deco,0];
    end
end
huff_code_ask = decode(mt_ask_deco,n,k,'hamming/binary');


fprintf("\n")
fprintf("接收方Hamming译码（ASK）\t\t:")
fprintf("%s", string(huff_code_ask))
fprintf("\n")

ask_decoded = huffmandeco(huff_code_ask, dict); % 根据huff_code和dict还原字符串
ask_decoded = string(ask_decoded);
fprintf("接收方Huffman译码（ASK）\t\t:")
fprintf("%s", string(ask_decoded))
fprintf("\n")

acc_ask = cal_acc(str, ask_decoded);
fprintf("ASK准确率\t\t\t\t\t:")
fprintf("%f", acc_ask)
fprintf("\n")


%------2PSK编码------------------------------------------------------------
mt_psk=[];
external_t=0+(T/sample):(T/sample):length(ham_code)+1;  % +1是为了防止相位+Pi后溢出
sin_wave = sin(2*pi/T .*external_t);
for i=1:length(ham_code)
    if(ham_code(i)==1)
        for j = 1:sample
            mt_psk = [mt_psk,sin_wave((i-1)*sample + j)];
        end
    else 
        for j = 1:sample
            mt_psk = [mt_psk,sin_wave((i-1)*sample + j + sample/T/2)];
        end
    end
end

mt_psk_noise = awgn(mt_psk, snr);  % 加噪声

if test_iter == 1
    subplot(4,1,3)
    plot(t,mt_psk);
    ylim([-1.5,1.5]);
    title('2PSK')
end

%------2PSK译码------------------------------------------------------------
carrier = sin(2*pi/T .*t);
% 课件上还有带通滤波器，带通实现起来不如低通容易，先暂时不加，看看效果
mt_psk_xx = carrier .* mt_psk;  % 与载波相乘后的psk信号
mt_psk_xx_filted = lowpass(mt_psk_xx,1/T/100,1/T);
mt_psk_deco = [];
for t_judge = 0.5:T:length(ham_code)-0.5 % t_sampling是抽样判决时刻
    if mt_psk_xx_filted(t_judge*sample)>0
        mt_psk_deco = [mt_psk_deco,1];
    else
        mt_psk_deco = [mt_psk_deco,0];
    end
end
huff_code_psk = decode(mt_psk_deco, n, k,'hamming/binary');


fprintf("\n")
fprintf("接收方Hamming译码（PSK）\t\t:")
fprintf("%s", string(huff_code_psk))
fprintf("\n")

psk_decoded = huffmandeco(huff_code_psk, dict); % 根据huff_code和dict还原字符串
psk_decoded = string(psk_decoded);
fprintf("接收方Huffman译码（PSK）\t\t:")
fprintf("%s", string(psk_decoded))
fprintf("\n")

acc_psk = cal_acc(str, psk_decoded);
fprintf("PSK准确率\t\t\t\t\t:")
fprintf("%f", acc_psk)
fprintf("\n")



%------2FSK编码------------------------------------------------------------
mt_fsk=[];
f_l = 1/T;
f_h = 2 * 1/T;
sin_wave_l = sin(2*pi*f_l .*t);
sin_wave_h = sin(2*pi*f_h .*t);
for i=1:length(ham_code)
    if(ham_code(i)==1)
        for j = 1:sample
            mt_fsk = [mt_fsk,sin_wave_l((i-1)*sample + j)];
        end
    else 
        for j = 1:sample
            mt_fsk = [mt_fsk,sin_wave_h((i-1)*sample + j)];
        end
    end
end

if test_iter == 1
    subplot(4,1,4)
    plot(t,mt_fsk);
    ylim([-1.5,1.5]);
    title('2FSK')
end

%------一些其他图像---------------------------------------------------------
if test_iter == 1
    figure(2)
    subplot(2,1,1)
    plot(t,mt_ask);
    ylim([-1.5,1.5]);
    title('2ASK')
    subplot(2,1,2)
    plot(t,mt_ask_noise);
    ylim([-1.5,1.5]);
    title('2ASK（噪声）')

    figure(3)
    subplot(2,1,1)
    plot(t,mt_psk);
    ylim([-1.5,1.5]);
    title('2PSK')
    subplot(2,1,2)
    plot(t,mt_psk_noise);
    ylim([-1.5,1.5]);
    title('2PSK（噪声）')
end

end


%------准确率判断-----------------------------------------------------------
function acc = cal_acc(x1, x2)
    acc = 0;
    if length(x1) ~= length(x2)
        fprintf("!!!字符长度不匹配!!!\n")
    else
        for i = 1:length(x1)
            if x1(i) == x2(i)
                acc = acc + 1/length(x1)*100;  % MATLAB没有 +=
            end
        end
    end
end
