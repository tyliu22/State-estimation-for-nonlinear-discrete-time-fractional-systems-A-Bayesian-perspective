%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order EKF
%   Ŀ�ģ�FCDKF��FEKF�����ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע�� ���棬����ʱ�估����
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clc
 clear
 tempp = 50;

 for kkk = 1:tempp

%���沽��
N = 50;

q = 1;                % ϵͳ������ֵ
r = 1;                % ����������ֵ
Q = 0.81;             % ϵͳ�����������
R = 0.25;             % ���������������

%GL�����¶̼���ԭ��ĳ���
L = N+1;

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%ϵͳ��������
I = eye(1,1);                %���ɵ�λ��

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%����
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

x_0  = 0;                    %��ʼ״̬     
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

X_esti = zeros(1,N);        %״̬���Ź���ֵ
P_xesti = zeros(1,N);       %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����

% ϵͳ�������������
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = f(X_real(1,k-1)) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = h(X_real(1,k)) + V_noise(1,k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------��������չ�������˲������ܲ���---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

for k=2:1:N
  %�������˲�
      %״̬Ԥ��:X_pre
        diff_X_esti = f(X_esti(1,k-1));
            %��������
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = 0;
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            end
        F = 6*cos(2*X_esti(1,k-1)) - exp(X_esti(1,k-1));
            
        P_xpred = (F-bino_fir(1,2))*P_xesti(1,k-1)*(F-bino_fir(1,2))'+ Q + rema_P;
        
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;
        
        %���㿨��������:Kk(2*1)
        H = 1;
        Kk = P_xpred*H'/(H*P_xpred*H' + R);
        
        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
        
        %���Ʒ���������
        P_xesti(1,k) = (I-Kk*H)*P_xpred;
end
FEKF_SISO_TIME(1,kkk)    = toc;
FEKF_ERROR_norm1(1,kkk)  = norm((X_real - X_esti),1);
FEKF_ERROR_norm2(1,kkk)  = norm((X_real - X_esti),2);
FEKF_ERROR_RMSE(1,kkk)   = sqrt(sum((X_real - X_esti).^2)/N);

 end



 FEKF_SISO_ERROR_TIME(1,:) = FEKF_ERROR_norm1(1,:);
 FEKF_SISO_ERROR_TIME(2,:) = FEKF_ERROR_norm2(1,:);
 FEKF_SISO_ERROR_TIME(3,:) = FEKF_SISO_TIME(1,:);
 FEKF_ERROR_norm1_average  = sum(FEKF_SISO_ERROR_TIME(1,:))/50
 FEKF_ERROR_norm2_average  = sum(FEKF_SISO_ERROR_TIME(2,:))/50
 FEKF_SISO_TIME_average    = sum(FEKF_SISO_ERROR_TIME(3,:))/50
%  save FEKF_SISO_ERROE_TIME1 FEKF_SISO_ERROR_TIME FEKF_ERROR_norm1_average ...
%       FEKF_ERROR_norm2_average FEKF_SISO_TIME_average





