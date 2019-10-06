%*************************************************************************%
%   �����������˲����渴��
%   ���ģ�     fractional order PF
%   Ŀ�ģ������������˲��㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й���
%
%   ��ע�������������˲����㷨����
%           ����ز���
%*************************************************************************%

clc
clear

% ָ����ǰ·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOSystemDataSet\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));
% SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));
Y_RealMeas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%���沽��
N = SystemParameter{1,1}{1};     % ���沽��
q = SystemParameter{1,1}{2};     % ϵͳ������ֵ
r = SystemParameter{1,1}{3};     % ����������ֵ
Q = SystemParameter{1,1}{4};     % ϵͳ��������
R = SystemParameter{1,1}{5};     % ������������
alpha = SystemParameter{1,1}{6}; % ϵͳ�״�

SimuTimes = N;
NumParticle = 10;  % ���Ӹ���

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);          %΢�ֽ״�Ϊ0.03ʱGL�����µ�ϵ��
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha(1,1)+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha(1,2)+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %΢�ֽ״�Ϊ0.1ʱGL�����µ�ϵ��
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha(1,3)+1)/(i-1))*bino_thi(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������
I = eye(3,3);                %���ɵ�λ��



X_RealState = zeros(3,SimuTimes);        % ϵͳ״̬��ʵֵ ��ʼֵ0

P_SampleCov = cell(1,SimuTimes);        % ��������
x_EstiState = zeros(3,SimuTimes);        % ״̬����ֵ
x_EstiState(:,1) = [0.1 0.1 0.1]';
P_SampleCov{1,1} = [10,0,0;0,10,0;0,0,10];   

% ��ʼ�����ֲ��ķ���

ParticleWeight    = zeros(1,NumParticle);     % ��ʼ��Ȩ��
x_SamplePart_temp = zeros(3,NumParticle);     % �м����
x_SampleParticle  = zeros(3,NumParticle,SimuTimes);

% Intinialization particle, prior distirbution p(x_0) 
for i = 1 : NumParticle
    % ��ʼ״̬���� x=0 ��ֵ������Ϊ sqrt(P) �ĸ�˹�ֲ�
    x_SampleParticle(:,i,1) = mvnrnd(x_EstiState(:,1) + q,P_SampleCov{1,1},1);
end


f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
   
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));


%% ����ʵ��״̬ calculate real state and measurement
for k = 2 : SimuTimes
    
 %% ����N������ 
 for i = 1 : NumParticle
     % Draw particle: x^i_k ~ p(x_k | x^i_k-1) state transform function
     
     % ������� Num_particle ������
     x_SamplePart_temp(:,i) = mvnrnd(f(x_SampleParticle(:,i,k-1)) + q,Q,1);
     % x_SamplePart_temp(1,i) =  f(x_SampleParticle(i,k-1)) + q + sqrt(Q) * randn;
     temp = [0;0;0];
         for j = 2 : 1 : k
            temp = temp + gamma{1,j}*x_SampleParticle(:,i,k+1-j);
         end
     x_SamplePart_temp(:,i) = x_SamplePart_temp(:,i) - temp;
     
     y_ParticleMeas = h(x_SamplePart_temp(:,i)) + r;     % ÿ�����Ӷ�Ӧ�Ĺ۲�ֵ
     ErrorMeas = Y_RealMeas(1,k) - y_ParticleMeas;   % ����ʵ�۲�֮�����Ȼ
     % Draw weight: w^i_k ~ p(z_k | x^i_k) measurement transform function
     % ����Ȩֵ������������й�
     %ParticleWeight(1,i) = h(x_SamplePart_temp(1,i)) + r + sqrt(R) * randn;
     ParticleWeight(1,i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-ErrorMeas^2 / 2 / R);
     % ÿ�����ӵ���Ȼ�����ƶ�
 end

%%
 % Ȩֵ��һ��
weight_sum = sum(ParticleWeight);
for i = 1 : NumParticle
    ParticleWeight(1,i) = ParticleWeight(1,i) / weight_sum;  % ��һ�����Ȩֵ q
end

 % ����Ȩֵ���²��� ����ز����㷨
 qtempsum = zeros(1,NumParticle);
  qtempsum(1,1) = ParticleWeight(1,1);
 for i = 2 : 1 : NumParticle
    qtempsum(1,i) = qtempsum(1,i-1) + ParticleWeight(1,i);
 end
 
  for i = 1 : NumParticle
      UniRandom = rand; % �������ȷֲ������
      for j = 1 : NumParticle
          % �ۼ�Ȩֵ
          %qtempsum = qtempsum + ParticleWeight(1,j);
          if qtempsum(1,j) >= UniRandom
              x_SampleParticle(:,i,k) = x_SamplePart_temp(:,j);
              break;
          %else
          %    x_SampleParticle(i,k) = x_SampleParticle(i,k-1);
          end
      end
  end
  
%% ����ϵͳ״̬

x_EstiState(:,k) = mean((x_SampleParticle(:,:,k))');

end

% LineWidth = 1.5;
% %%
% t = 1 : SimuTimes;
% figure;
% plot(t, X_RealState, 'r', t, x_EstiState, 'b--','linewidth',LineWidth);
% Esitimated_state = legend('Real Value','Estimated Value','Location','best');
% set(Esitimated_state,'Interpreter','latex')
% set(gcf,'Position',[200 200 400 300]); 
% axis([0 50 -6 6]) % ������������ָ��������
% axis normal
% set(gca,'FontSize',10); 
% xlabel('time step','FontSize',7); 
% ylabel('state','FontSize',7);
% % ����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)

FPF_X_esti_100 = x_EstiState;

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOEstimatedState\';
save(strcat(path,'FPF_EstimatedState_100','.mat'),'FPF_X_esti_100');



