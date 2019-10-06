%*************************************************************************%
%                  分数阶无迹卡尔曼滤波器仿真复现                         %
%_________________________________________________________________________%
%   论文 : fractional order FUKF
%   目的 : 分数阶无迹卡尔曼滤波器仿真复现
%   函数实验 :
%               D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                       y_k = x_k + v_k
%
%   结果 : 
%
%   备注 : 
%
%*************************************************************************%

% save和load的使用方法：
% save(strcat('D:\存储文件夹名\',filesep,'example','.mat'),'data');
% load(strcat('D:\存储文件夹名\',filesep,'example','.mat'));

% SystemParameter = zeros(1,6);
% SystemParameter(1) = N;     % 仿真步长
% SystemParameter(2) = q;     % 系统噪声均值
% SystemParameter(3) = r;     % 测量噪声均值
% SystemParameter(4) = Q;     % 系统噪声方差
% SystemParameter(5) = R;     % 测量噪声方差
% SystemParameter(6) = alpha; % 系统阶次

clc
clear

% 指定当前路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';
SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));

%仿真步长
N = SystemParameter(1);     % 仿真步长
q = SystemParameter(2);     % 系统噪声均值
r = SystemParameter(3);     % 测量噪声均值
Q = SystemParameter(4);     % 系统噪声方差
R = SystemParameter(5);     % 测量噪声方差
alpha = SystemParameter(6); % 系统阶次

Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

 % GL 定义下短记忆原理的长度
 L = N+1;

% 计算 alpha 阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       % 微分阶次为 alpha 时 GL 定义下的系数
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

I = eye(1,1);                 % 生成单位阵

X_esti = zeros(1,N);          % 状态最优估计值

% 系统函数与测量函数
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

%*************************************************************************%
%-------------------------无迹卡尔曼滤波器性能测试------------------------%
%*************************************************************************%

P_xesti = cell(1,N);            % 估计误差方差阵

% 初始值设置（初始矩阵不能为零）
P_pred_0     = 10;              % 初始预测方差阵
P_xesti{1,1} = P_pred_0;        % 初始估计方差阵

state_dim = 1;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(1, L_sample);
SigmaWeight = zeros(1, L_sample);
%GammaPoints = zeros(1, L_sample);
ChiPoints   = zeros(1, L_sample);

for k=2:1:N
    % 对称采样
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_esti(:,k-1), P_xesti{1,k-1}); 
    
    for i = 1 : 1 : L_sample
        ChiPoints(:,i) = f(SigmaPoints(:,i)) + q;
    end

    % Predicted state
        X_pre = 0;
    for i = 1 : 1 : L_sample
        X_pre = X_pre +  SigmaWeight(:,i)*ChiPoints(:,i);
    end
        % 计算余项
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
    X_pre = X_pre - rema;

    % Predicted state error covariance 
    P_xpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_xpre = P_xpre + ...
                 SigmaWeight(:,i)*(ChiPoints(:,i)-X_pre)*(ChiPoints(:,i)-X_pre)' + ...
                 SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(1,k-1))*(ChiPoints(:,i)-X_pre)'*bino_fir(1,2)' ...
                 + bino_fir(1,2)*SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(1,k-1))*(ChiPoints(:,i)-X_pre)';
    end
        % 计算余项
        rema_P = 0;
        if k>L+1
            for i = 2:1:L+2
                rema_P = rema_P + bino_fir(1,i)*P_xesti{1,k+1-i}*bino_fir(1,i)';
            end
        else
            for i = 2:1:k
                rema_P = rema_P + bino_fir(1,i)*P_xesti{1,k+1-i}*bino_fir(1,i)';
            end
        end
    P_xpre = P_xpre + rema_P + Q;

    % measurement update
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_pre, P_xpre); 
    
    for i = 1 : 1 : L_sample
        ChiPoints(:,i) = h(SigmaPoints(:,i)) + r;
    end
    
    % Predicted measurement
    Z_pre = 0;
    for i = 1 : 1 : L_sample
        Z_pre = Z_pre + SigmaWeight(:,i)*ChiPoints(:,i);
    end
    
    % Predicted measurement error covariance 
    P_zpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_zpre = P_zpre + SigmaWeight(:,i)*(ChiPoints(:,i)-Z_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    P_zpre = P_zpre + R;
   
    % cross-variance 
    P_xzpre = 0;
    for i = 1 : 1 : L_sample 
        P_xzpre = P_xzpre + SigmaWeight(:,i)*(SigmaPoints(:,i)-X_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    
    % Kalman gain
    Kk = P_xzpre/P_zpre;
    
    % estimated state
    X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_pre );

    % estimation error covariance
    P_xesti{1,k} = P_xpre - Kk*P_zpre*Kk';

end


% % 输入与测量输出图
% k = 1:1:N;
% LineWidth = 1.5;
% 
% % square error
% figure;
% plot(k,X_state_real(1,:),'r',k,X_state_esti(1,:),'b--','linewidth',LineWidth);
% % set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
%  axis normal
%  axis([ 0 N -6 6 ])
% ylabel('x','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% % 设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% legend('Real state 1','Estimation state 1','Location','best');


FUKF_X_esti = X_esti;

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
save(strcat(path,'FUKF_EstimatedState','.mat'),'FUKF_X_esti');


