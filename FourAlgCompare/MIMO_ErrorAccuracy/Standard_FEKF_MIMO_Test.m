%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order EKF
%   目的：FCDKF与FEKF的性能比较
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：
%
%   备注：
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOSystemDataSet\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));
% SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));
Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%仿真步长
N = SystemParameter{1,1}{1};     % 仿真步长
q = SystemParameter{1,1}{2};     % 系统噪声均值
r = SystemParameter{1,1}{3};     % 测量噪声均值
Q = SystemParameter{1,1}{4};     % 系统噪声方差
R = SystemParameter{1,1}{5};     % 测量噪声方差
alpha = SystemParameter{1,1}{6}; % 系统阶次

%GL定义下短记忆原理的长度
L = N+1;

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);          %微分阶次为0.03时GL定义下的系数
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha(1,1)+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %微分阶次为1.2时GL定义下的系数
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha(1,2)+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %微分阶次为0.1时GL定义下的系数
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha(1,3)+1)/(i-1))*bino_thi(1,i-1);  
end

%计算GL微分定义下系数矩阵
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%系统矩阵设置
I = eye(3,3);                %生成单位阵

% 系统函数与测量函数
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------分数阶扩展卡尔曼滤波器性能测试---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti = zeros(3,N);          %状态最优估计值
X_esti(:,1) = [0.1 0.1 0.1]'; %初始状态估计
P_xesti = cell(1,N);         %估计误差方差阵

%初始值设置（初始矩阵不能为零）
P_pred_0 = [10,0,0;0,10,0;0,0,10];   %初始预测方差阵
P_xesti{1,1} = P_pred_0;             %初始估计方差阵

for k=2:1:N
  %卡尔曼滤波
      %状态预测:X_pre
        diff_X_esti = f(X_esti(:,k-1));
            %计算余项
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
            rema_P = [0,0,0;0,0,0;0,0,0];
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_xesti{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + gamma{1,i}*P_xesti{1,k+1-i}*gamma{1,i}';
                end
            end
        F = [0.2*X_esti(3,k-1), -0.3*sin(X_esti(2,k-1)), 0.2*X_esti(1,k-1); ...
            -0.2*sin(X_esti(1,k-1)), 0, - 0.5*cos(X_esti(3,k-1)); ...
            0, 0.7*cos(X_esti(2,k-1)) * cos(X_esti(3,k-1)), -0.7*sin(X_esti(2,k-1))*sin(X_esti(3,k-1))];        
            
        P_xpred = (F-gamma{1,2})*P_xesti{1,k-1}*(F-gamma{1,2})'+ Q + rema_P;
        
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;
        
        %计算卡尔曼增益:Kk(2*1)
        H = [0.3*X_pre(2), 0.3*X_pre(1), - 0.2*cos(X_pre(3)) ];    

        Kk = P_xpred*H'/(H*P_xpred*H' + R);
        
        %状态更新
        X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
        
        %估计方差矩阵更新
        P_xesti{1,k} = (I-Kk*H)*P_xpred;
end

% %输入与测量输出图
% k = 1:1:N;
% 
% LineWidth = 1.5;
% 
% figure;
% plot(k,X_real(1,:),'r',k,X_esti(1,:),'b--','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
% axis normal
% axis([0 N -6 6 ])
% ylabel('x','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% % 设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');

FEKF_X_esti = X_esti;

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOEstimatedState\';
save(strcat(path,'FEKF_EstimatedState','.mat'),'FEKF_X_esti');
