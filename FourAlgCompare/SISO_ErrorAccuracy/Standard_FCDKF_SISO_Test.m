%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order CDKF
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

%GL定义下短记忆原理的长度
L = N+1;
h_con = sqrt(3);

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%系统矩阵设置
I = eye(1,1);                %生成单位阵

% 系统函数与测量函数
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------分数阶卡尔曼滤波器性能测试------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti  = zeros(1,N);        %状态最优估计值
P_xesti = zeros(1,N);        %估计误差方差阵

P_pred_0 = 20;              %初始预测方差阵
P_xesti(1,1) = P_pred_0;     %初始估计方差阵

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

            %对误差矩阵进行cholsky分解
            S_chol = chol(P_xesti(1,k-1))';

            %计算余项
            rema_P = 0;
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            end

        %临时变量 temp_fun : 函数差值,函数为单变量
        temp_fun = f(X_esti(1,k-1)+h_con*S_chol) - f(X_esti(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q;
        
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;

        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %状态更新
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %估计方差矩阵更新
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';
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

FCDKF_X_esti = X_esti;

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
save(strcat(path,'FCDKF_EstimatedState','.mat'),'FCDKF_X_esti');





