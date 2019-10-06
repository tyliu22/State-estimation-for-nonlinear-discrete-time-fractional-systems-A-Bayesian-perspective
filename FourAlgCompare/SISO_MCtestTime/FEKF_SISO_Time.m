%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order EKF
%   目的：FCDKF与FEKF的性能比较
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                              y_k = x_k + v_k
%   结果：
%
%   备注： 仿真，计算时间及误差范数
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clc
 clear
 tempp = 50;

 for kkk = 1:tempp

%仿真步长
N = 50;

q = 1;                % 系统噪声均值
r = 1;                % 测量噪声均值
Q = 0.81;             % 系统噪声方差矩阵
R = 0.25;             % 测量噪声方差矩阵

%GL定义下短记忆原理的长度
L = N+1;

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%系统矩阵设置
I = eye(1,1);                %生成单位阵

%状态测量初始化
X_real = zeros(1,N);         %真实状态
Z_meas = zeros(1,N);         %实际观测值

%噪声
W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

x_0  = 0;                    %初始状态     
X_real(1,1) = x_0;           %真实状态初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

X_esti = zeros(1,N);        %状态最优估计值
P_xesti = zeros(1,N);       %估计误差方差阵

%初始值设置（初始矩阵不能为零）
P_pred_0 = 100;              %初始预测方差阵
P_xesti(1,1) = P_pred_0;     %初始估计方差阵

% 系统函数与测量函数
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

for k=2:1:N
    %计算实际状态
    diff_X_real = f(X_real(1,k-1)) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %实际观测值
    Z_meas(1,k) = h(X_real(1,k)) + V_noise(1,k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------分数阶扩展卡尔曼滤波器性能测试---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

for k=2:1:N
  %卡尔曼滤波
      %状态预测:X_pre
        diff_X_esti = f(X_esti(1,k-1));
            %计算余项
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
        X_pre = diff_X_esti - rema + q;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
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
        
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;
        
        %计算卡尔曼增益:Kk(2*1)
        H = 1;
        Kk = P_xpred*H'/(H*P_xpred*H' + R);
        
        %状态更新
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
        
        %估计方差矩阵更新
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





