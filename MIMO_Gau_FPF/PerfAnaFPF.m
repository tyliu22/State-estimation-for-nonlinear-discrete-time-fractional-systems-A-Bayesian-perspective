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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_Gau_FPF\MIMOEstimatedState\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));

%仿真步长
N = SystemParameter{1,1}{1};     % 仿真步长
q = SystemParameter{1,1}{2};     % 系统噪声均值
r = SystemParameter{1,1}{3};     % 测量噪声均值
Q = SystemParameter{1,1}{4};     % 系统噪声方差
R = SystemParameter{1,1}{5};     % 测量噪声方差
alpha = SystemParameter{1,1}{6}; % 系统阶次


FPF_X_esti_500 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_500','.mat'))));
FPF_X_esti_250 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_250','.mat'))));
FPF_X_esti_100 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_100','.mat'))));

RealState = cell2mat(struct2cell(load(strcat(path,'RealState','.mat'))));


% FEKF_error = FEKF_X_esti - RealState;


k = 1:1:50;
LineWidth = 1.5;
FontSize = 7;

%% error
figure;
subplot(311)
plot(k,RealState(1,:),'b-',k,FPF_X_esti_250(1,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 -0.2 0.5])
ylabel('$x_1$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h1 = legend('$x_1$','$\hat{x}_1$','Location','best');
set(h1,'Interpreter','latex')

subplot(312)
plot(k,RealState(2,:),'b-',k,FPF_X_esti_250(2,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 -1 1.2])
ylabel('$x_2$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h2 = legend('$x_2$','$\hat{x}_2$','Location','best');
set(h2,'Interpreter','latex')

subplot(313)
plot(k,RealState(3,:),'b-',k,FPF_X_esti_250(3,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 -0.5 0.8])
ylabel('$x_3$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h3 = legend('$x_3$','$\hat{x}_3$','Location','best');
set(h3,'Interpreter','latex')



%% RMSE
% SE
Error_FPF_X_esti_100 = (FPF_X_esti_100(:,:) - RealState(:,:)).^2; % .^2
Error_FPF_X_esti_250 = (FPF_X_esti_250(:,:) - RealState(:,:)).^2;
Error_FPF_X_esti_500 = (FPF_X_esti_500(:,:) - RealState(:,:)).^2;


%% error
figure;
subplot(311)
plot(k,Error_FPF_X_esti_100(1,:),'b-',k,Error_FPF_X_esti_250(1,:),'r--',k,Error_FPF_X_esti_500(1,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 0 0.3])
ylabel('SE($x_1$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he1 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he1,'Interpreter','latex')

subplot(312)
plot(k,Error_FPF_X_esti_100(2,:),'b-',k,Error_FPF_X_esti_250(2,:),'r--',k,Error_FPF_X_esti_500(2,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 0 4])
ylabel('SE($x_2$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he2 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he2,'Interpreter','latex')

subplot(313)
plot(k,Error_FPF_X_esti_100(3,:),'b-',k,Error_FPF_X_esti_250(3,:),'r--',k,Error_FPF_X_esti_500(3,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  设置坐标轴在指定的区间
axis normal
axis([0 25 0 0.5])
ylabel('SE($x_3$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% 设置坐标轴刻度字体名称，大小
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he3 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he3,'Interpreter','latex')

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

