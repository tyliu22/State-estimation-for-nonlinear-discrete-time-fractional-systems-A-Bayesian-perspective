%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文： 
%   目的：FCDKF, FEKF, FPF, FUKF 估计精度性能比较
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

% SystemParameter(1) = N;     % 仿真步长
% SystemParameter(2) = q;     % 系统噪声均值
% SystemParameter(3) = r;     % 测量噪声均值
% SystemParameter(4) = Q;     % 系统噪声方差
% SystemParameter(5) = R;     % 测量噪声方差
% SystemParameter(6) = alpha; % 系统阶次

% FPF_EstimatedState   FPF_X_esti 
% FUKF_EstimatedState  FUKF_X_esti 
% FEKF_EstimatedState  FEKF_X_esti
% FCDKF_EstimatedState FCDKF_X_esti

clc
clear

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOEstimatedState\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));

%仿真步长
N = SystemParameter{1,1}{1};     % 仿真步长
q = SystemParameter{1,1}{2};     % 系统噪声均值
r = SystemParameter{1,1}{3};     % 测量噪声均值
Q = SystemParameter{1,1}{4};     % 系统噪声方差
R = SystemParameter{1,1}{5};     % 测量噪声方差
alpha = SystemParameter{1,1}{6}; % 系统阶次


FPF_X_esti = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState','.mat'))));

RealState = cell2mat(struct2cell(load(strcat(path,'RealState','.mat'))));


% FEKF_error = FEKF_X_esti - RealState;


k = 1:1:50;
LineWidth = 1.5;
FontSize = 8;

%% error
figure;
subplot(311)
plot(k,RealState(1,:),'b-',k,FPF_X_esti(1,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 0 0.5])
ylabel('$x_1$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h1 = legend('$x_1$','$\hat{x}_1$','Location','best');
set(h1,'Interpreter','latex')

subplot(312)
plot(k,RealState(2,:),'b-',k,FPF_X_esti(2,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 -0.5 1])
ylabel('$x_2$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h2 = legend('$x_2$','$\hat{x}_2$','Location','best');
set(h2,'Interpreter','latex')

subplot(313)
plot(k,RealState(3,:),'b-',k,FPF_X_esti(3,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 -0.2 0.8])
ylabel('$x_3$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h3 = legend('$x_3$','$\hat{x}_3$','Location','best');
set(h3,'Interpreter','latex')

%% RMSE
% SE
Error_FPF_X_esti   = (FPF_X_esti   - RealState).^2;
Error_FUKF_X_esti  = (FUKF_X_esti  - RealState).^2;
Error_FEKF_X_esti  = (FEKF_X_esti  - RealState).^2;
Error_FCDKF_X_esti = (FCDKF_X_esti - RealState).^2;

RMSE_FPF_X_esti   = zeros(1,N);
RMSE_FUKF_X_esti  = zeros(1,N);
RMSE_FEKF_X_esti  = zeros(1,N);
RMSE_FCDKF_X_esti = zeros(1,N);

RMSE_FPF_X_esti(1,1)   = Error_FPF_X_esti(1,1);
RMSE_FUKF_X_esti(1,1)  = Error_FUKF_X_esti(1,1);
RMSE_FEKF_X_esti(1,1)  = Error_FEKF_X_esti(1,1);
RMSE_FCDKF_X_esti(1,1) = Error_FCDKF_X_esti(1,1);

 for i = 2:1:N
     RMSE_FPF_X_esti(1,i)   = RMSE_FPF_X_esti(1,i-1)   + Error_FPF_X_esti(1,i);
     RMSE_FUKF_X_esti(1,i)  = RMSE_FUKF_X_esti(1,i-1)  + Error_FUKF_X_esti(1,i);
     RMSE_FEKF_X_esti(1,i)  = RMSE_FEKF_X_esti(1,i-1)  + Error_FEKF_X_esti(1,i);
     RMSE_FCDKF_X_esti(1,i) = RMSE_FCDKF_X_esti(1,i-1) + Error_FCDKF_X_esti(1,i);
 end
 
  for i = 1:1:N
     RMSE_FPF_X_esti(1,i)   = sqrt( RMSE_FPF_X_esti(1,i) / i );
     RMSE_FUKF_X_esti(1,i)  = sqrt( RMSE_FUKF_X_esti(1,i) / i );
     RMSE_FEKF_X_esti(1,i)  = sqrt( RMSE_FEKF_X_esti(1,i) / i );
     RMSE_FCDKF_X_esti(1,i) = sqrt( RMSE_FCDKF_X_esti(1,i) / i );
  end
 
figure;
plot(k,RMSE_FEKF_X_esti(1,:),'r-',k,RMSE_FCDKF_X_esti(1,:),'b--', ...
     k,RMSE_FUKF_X_esti(1,:),'g:', k,RMSE_FPF_X_esti(1,:),'c-.', 'linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 N 0 1.5])
ylabel('RMSE','FontSize',FontSize)
xlabel('iteration times','FontSize',FontSize)
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',FontSize)
legend('FEKF','FCDKF','FUKF','FPF','Location','best');
 
 
 
 





%% Square Error
figure;
plot(k,Error_FEKF_X_esti(1,:),'r-',k,Error_FCDKF_X_esti(1,:),'b--', ...
     k,Error_FUKF_X_esti(1,:),'g:', k,Error_FPF_X_esti(1,:),'c-.', 'linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 N 0 7])
ylabel('SE','FontSize',FontSize)
xlabel('iteration times','FontSize',FontSize)
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',FontSize)
legend('FEKF','FCDKF','FUKF','FPF','Location','best');


%% 
% 'r' 红色 'm' 粉红
% 'g' 绿色 'c' 青色
% 'b' 兰色 'w' 白色
% 'y' 黄色 'k' 黑色

% '-' 实线 '--' 虚线
% ':' 点线 '-.' 点划线

% '.' 用点号绘制各数据点 '^' 用上三角绘制各数据点
% '+' 用'+'号绘制各数据点 'v' 用下三角绘制各数据点
% '*' 用'*'号绘制各数据点 '>' 用右三角绘制各数据点
% ' .' 用'.'号绘制各数据点 '<' 用左三角绘制各数据点
% 's'或squar 用正方形绘制各数据点'p' 用五角星绘制各数据点
% 'd'或diamond用菱形绘制各数据点 'h' 用六角星绘制各数据点

