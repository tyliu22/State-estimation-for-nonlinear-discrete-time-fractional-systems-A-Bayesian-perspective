%*************************************************************************%
%   分数阶粒子滤波仿真复现
%   论文：     fractional order PF
%   目的：分数阶粒子滤波算法测试
%         对系统噪声均值进行估计
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：较好的对状态进行估计
%
%   备注：分数阶粒子滤波的算法测试
%           随机重采样
%*************************************************************************%

clc;
clear all;

% 指定当前路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';
SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));

%仿真步长
SimuTimes = SystemParameter(1);     % 仿真步长
q = SystemParameter(2);     % 系统噪声均值
r = SystemParameter(3);     % 测量噪声均值
Q = SystemParameter(4);     % 系统噪声方差
R = SystemParameter(5);     % 测量噪声方差
alpha = SystemParameter(6); % 系统阶次

Y_RealMeas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

NumParticle = 500;  % 粒子个数

%系统矩阵设置
I = eye(1,1);               % 生成单位阵

X_RealState = zeros(1,SimuTimes); % 系统状态真实值 初始值0

P_SampleCov = zeros(1,SimuTimes);        % 采样方差
x_EstiState = zeros(1,SimuTimes);        % 状态估计值
P_SampleCov(1,1) = 20;                   % 初始采样分布的方差

ParticleWeight    = zeros(1,NumParticle);     % 初始化权重
x_SamplePart_temp = zeros(1,NumParticle);     % 中间变量
x_SampleParticle  = zeros(NumParticle,SimuTimes);

% Intinialization particle, prior distirbution p(x_0) 
for i = 1 : NumParticle
    % 初始状态服从 x=0 均值，方差为 sqrt(P) 的高斯分布
    x_SampleParticle(i,1) = x_EstiState(1,1) + q + sqrt(P_SampleCov(1,1)) * randn; 
end

f = @(x)3*sin(2*x) - exp(x);
h = @(x)x;

% 计算alpha阶次对应的GL定义系数 binomial coefficient
bino_fir = zeros(1,SimuTimes);       % 微分阶次为0.7时GL定义下的系数
bino_fir(1,1) = 1;
for i = 2:1:NumParticle
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end


%% 计算实际状态 calculate real state and measurement
for k = 2 : SimuTimes
    
 %% 采样N个粒子 
 for i = 1 : NumParticle
     % Draw particle: x^i_k ~ p(x_k | x^i_k-1) state transform function
     % 采样获得 Num_particle 个粒子
     x_SamplePart_temp(1,i) =  f(x_SampleParticle(i,k-1)) + q + sqrt(Q) * randn;
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + bino_fir(1,j)*x_SampleParticle(i,k+1-j);
         end
     x_SamplePart_temp(1,i) = x_SamplePart_temp(1,i) - temp;
     y_ParticleMeas = h(x_SamplePart_temp(1,i)) + r;     % 每个粒子对应的观测值
     ErrorMeas = Y_RealMeas(1,k) - y_ParticleMeas;   % 与真实观测之间的似然
     % Draw weight: w^i_k ~ p(z_k | x^i_k) measurement transform function
     % 粒子权值，与测量方程有关
     %ParticleWeight(1,i) = h(x_SamplePart_temp(1,i)) + r + sqrt(R) * randn;
     ParticleWeight(1,i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-ErrorMeas^2 / 2 / R);
     % 每个粒子的似然即相似度
 end

%%
 % 权值归一化
weight_sum = sum(ParticleWeight);
for i = 1 : NumParticle
    ParticleWeight(1,i) = ParticleWeight(1,i) / weight_sum;  % 归一化后的权值 q
end

 % 根据权值重新采样 随机重采样算法
 qtempsum = zeros(1,NumParticle);
  qtempsum(1,1) = ParticleWeight(1,1);
 for i = 2 : 1 : NumParticle
    qtempsum(1,i) = qtempsum(1,i-1) + ParticleWeight(1,i);
 end
 
  for i = 1 : NumParticle
      UniRandom = rand; % 产生均匀分布随机数
      for j = 1 : NumParticle
          % 累计权值
          %qtempsum = qtempsum + ParticleWeight(1,j);
          if qtempsum(1,j) >= UniRandom
              x_SampleParticle(i,k) = x_SamplePart_temp(1,j);
              break;
          %else
          %    x_SampleParticle(i,k) = x_SampleParticle(i,k-1);
          end
      end
  end
  
%% 估计系统状态

x_EstiState(1,k) = mean(x_SampleParticle(:,k));

end

% LineWidth = 1.5;
% %%
% t = 1 : SimuTimes;
% figure;
% plot(t, X_RealState, 'r', t, x_EstiState, 'b--','linewidth',LineWidth);
% Esitimated_state = legend('Real Value','Estimated Value','Location','best');
% set(Esitimated_state,'Interpreter','latex')
% set(gcf,'Position',[200 200 400 300]); 
% axis([0 50 -6 6]) % 设置坐标轴在指定的区间
% axis normal
% set(gca,'FontSize',10); 
% xlabel('time step','FontSize',7); 
% ylabel('state','FontSize',7);
% % 设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)

FPF_X_esti = x_EstiState;

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
save(strcat(path,'FPF_EstimatedState','.mat'),'FPF_X_esti');



