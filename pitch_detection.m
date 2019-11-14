clc;
clear all;
close all;        

% 读取语料
[y,fs,wmode,fidx]=readwav('number.wav','p',-1,-1);
figure(1);stem(y,'.');title('原始音频波形');%显示声音的波形
%l=length(y);

%分帧
%f=enframe(y,hamming(160));                   %分帧长:20ms(160样点)
%figure(2);stem(f,'.');title('分帧后的音频波形');%显示声音的波形
nFrames = floor(length(y)/160) ; % 总帧数

%加窗，计算短时能量
w=window(@hamming,160);
for k = 1:nFrames
    idy = (k-1) * 160+(1:160);
    y_sub = y(idy) .* w;
    E(k) = sum(y_sub.^2); 
end
figure(2);stem(E,'.');title('全区间的短时能量图');

%如何设定能量阈值为0.001得到安静、有声的区间？
for i = 1:nFrames
    if E(i)<0.001
        E(i)=0
    end
end
figure(3);stem(E,'.');title('有声区间的短时能量图');%处理后的帧的短时能量图像

%计算短时自相关
n=160;
for m=1:length(y)/n            %对每一帧求短时自相关函数，每帧的Rm最大值存在N(m)里
    for k=20:100
        Rm(k)=0;
        for i=1:160
            Rm(k)=Rm(k)+y(i+(m-1)*n)*y(i+k+(m-1)*n);
        end
    end
    %p=Rm(10:100);                %防止误判，去掉前边10个数值较大的点
    p=Rm;
    [Rmax,N(m)]=max(p);        %读取第一个自相关函数的最大点  找到A中那些最大值的索引位置，将他们放在向量N中返回。
end                            

%计算基音周期和基音频率
%N=N+10;                        %补回前边去掉的10个点
%T=N/8;                         %算出对应的周期
T=(1/8)*N                       %f=8kHz 对应周期为（1/f）*N 对应单位为ms
figure(4);stem(T,'.');axis([0 length(T) 0 20]);
xlabel('帧数(n)');ylabel('周期(ms)');title('初始未处理的基音周期');

T1= medfilt1(T,7);             %去除野点，中值平滑
for k = 1:nFrames
    if E(k)==0
        T1(k)=0;
    end
end
figure(5);stem(T1,'.');axis([0 length(T1) 0 20]);
xlabel('帧数(n)');ylabel('周期(ms)');title('中值滤波后的基音周期');

F1= 1000./T1;                  %计算基音频率（对于能量为0的无声区间设置基音频率值为0）
for k = 1:nFrames
    if E(k)==0
        F1(k)=0;
    end
end
figure(6);stem(F1,'.');
xlabel('帧数(n)');ylabel('频率(Hz)');title('基音频率');

%计算数字“0”的差分方程系数和预测增益
p=10;
w=window(@hamming,160);  
y1=y(99*160:100*160-1);  %取数字“0”的一帧（第100帧），进行汉明加窗
A=lpc(y1.*w,p);  %得到系数（11个系数中第一个系数值为1）
est_Frame=filter([0 -A(2:end)],1,y1);%estimate frame(lp)预测第100/101帧?
FFT_est=fft(est_Frame);
%y2=y(100*160:101*160-1);
err=y1-est_Frame;%calculate the residual signal.
prodictive_gain=y1/err;
figure(7);
%subplot(221);plot(1:3360,y1,1:3360,est_Frame,'-r');grid;title('原始语音帧 vs.预测后的语音帧');
subplot(221);plot(y1);grid;title('原始语音段');
%subplot(222);plot(est_Frame);grid;title('预测语音段');
subplot(222);plot(est_Frame);grid;title({'预测语音段及10个预测系数',[num2str(A(2)),',',num2str(A(3)),',',num2str(A(4)),',',num2str(A(5)),',',num2str(A(6)),',',num2str(A(7)),',',num2str(A(8)),',',num2str(A(9)),',',num2str(A(10)),',',num2str(A(11))]});
subplot(223);plot(err);grid;title('误差');
subplot(224);plot(prodictive_gain);grid;title('预测增益');
