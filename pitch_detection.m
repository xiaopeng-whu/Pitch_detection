clc;
clear all;
close all;        

% ��ȡ����
[y,fs,wmode,fidx]=readwav('number.wav','p',-1,-1);
figure(1);stem(y,'.');title('ԭʼ��Ƶ����');%��ʾ�����Ĳ���
%l=length(y);

%��֡
%f=enframe(y,hamming(160));                   %��֡��:20ms(160����)
%figure(2);stem(f,'.');title('��֡�����Ƶ����');%��ʾ�����Ĳ���
nFrames = floor(length(y)/160) ; % ��֡��

%�Ӵ��������ʱ����
w=window(@hamming,160);
for k = 1:nFrames
    idy = (k-1) * 160+(1:160);
    y_sub = y(idy) .* w;
    E(k) = sum(y_sub.^2); 
end
figure(2);stem(E,'.');title('ȫ����Ķ�ʱ����ͼ');

%����趨������ֵΪ0.001�õ����������������䣿
for i = 1:nFrames
    if E(i)<0.001
        E(i)=0
    end
end
figure(3);stem(E,'.');title('��������Ķ�ʱ����ͼ');%������֡�Ķ�ʱ����ͼ��

%�����ʱ�����
n=160;
for m=1:length(y)/n            %��ÿһ֡���ʱ����غ�����ÿ֡��Rm���ֵ����N(m)��
    for k=20:100
        Rm(k)=0;
        for i=1:160
            Rm(k)=Rm(k)+y(i+(m-1)*n)*y(i+k+(m-1)*n);
        end
    end
    %p=Rm(10:100);                %��ֹ���У�ȥ��ǰ��10����ֵ�ϴ�ĵ�
    p=Rm;
    [Rmax,N(m)]=max(p);        %��ȡ��һ������غ���������  �ҵ�A����Щ���ֵ������λ�ã������Ƿ�������N�з��ء�
end                            

%����������ںͻ���Ƶ��
%N=N+10;                        %����ǰ��ȥ����10����
%T=N/8;                         %�����Ӧ������
T=(1/8)*N                       %f=8kHz ��Ӧ����Ϊ��1/f��*N ��Ӧ��λΪms
figure(4);stem(T,'.');axis([0 length(T) 0 20]);
xlabel('֡��(n)');ylabel('����(ms)');title('��ʼδ����Ļ�������');

T1= medfilt1(T,7);             %ȥ��Ұ�㣬��ֵƽ��
for k = 1:nFrames
    if E(k)==0
        T1(k)=0;
    end
end
figure(5);stem(T1,'.');axis([0 length(T1) 0 20]);
xlabel('֡��(n)');ylabel('����(ms)');title('��ֵ�˲���Ļ�������');

F1= 1000./T1;                  %�������Ƶ�ʣ���������Ϊ0�������������û���Ƶ��ֵΪ0��
for k = 1:nFrames
    if E(k)==0
        F1(k)=0;
    end
end
figure(6);stem(F1,'.');
xlabel('֡��(n)');ylabel('Ƶ��(Hz)');title('����Ƶ��');

%�������֡�0���Ĳ�ַ���ϵ����Ԥ������
p=10;
w=window(@hamming,160);  
y1=y(99*160:100*160-1);  %ȡ���֡�0����һ֡����100֡�������к����Ӵ�
A=lpc(y1.*w,p);  %�õ�ϵ����11��ϵ���е�һ��ϵ��ֵΪ1��
est_Frame=filter([0 -A(2:end)],1,y1);%estimate frame(lp)Ԥ���100/101֡?
FFT_est=fft(est_Frame);
%y2=y(100*160:101*160-1);
err=y1-est_Frame;%calculate the residual signal.
prodictive_gain=y1/err;
figure(7);
%subplot(221);plot(1:3360,y1,1:3360,est_Frame,'-r');grid;title('ԭʼ����֡ vs.Ԥ��������֡');
subplot(221);plot(y1);grid;title('ԭʼ������');
%subplot(222);plot(est_Frame);grid;title('Ԥ��������');
subplot(222);plot(est_Frame);grid;title({'Ԥ�������μ�10��Ԥ��ϵ��',[num2str(A(2)),',',num2str(A(3)),',',num2str(A(4)),',',num2str(A(5)),',',num2str(A(6)),',',num2str(A(7)),',',num2str(A(8)),',',num2str(A(9)),',',num2str(A(10)),',',num2str(A(11))]});
subplot(223);plot(err);grid;title('���');
subplot(224);plot(prodictive_gain);grid;title('Ԥ������');
