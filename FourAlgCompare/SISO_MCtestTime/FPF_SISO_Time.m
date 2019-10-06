%*************************************************************************%
%   �����������˲����渴��
%   ���ģ�     fractional order PF
%   Ŀ�ģ������������˲��㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й���
%
%   ��ע�� ���棬����ʱ�估����
%           
%*************************************************************************%

clc
clear
tempp = 50;

 for kkk = 1:tempp

%���沽��
SimuTimes = 50;    % ����ʱ��

LineWidth = 1.5;
NumParticle = 100;  % ���Ӹ���

%ϵͳ��������
I = eye(1,1);               % ���ɵ�λ��

%����
q = 1;                % ϵͳ������ֵ
r = 1;                % ����������ֵ
Q = 0.81;             % ϵͳ�����������
R = 0.25;             % ���������������

W_noise = sqrt(Q)*randn(1,SimuTimes) + q;    % ϵͳ����
V_noise = sqrt(R)*randn(1,SimuTimes) + r;    % ��������

X_RealState = zeros(1,SimuTimes); % ϵͳ״̬��ʵֵ ��ʼֵ0
Y_RealMeas = zeros(1,SimuTimes);  % ϵͳ״̬��ʵֵ ��ʼֵ0
Y_RealMeas(1,1) = X_RealState(1,1) + sqrt(R) * randn;

P_SampleCov = zeros(1,SimuTimes);        % ��������
x_EstiState = zeros(1,SimuTimes);        % ״̬����ֵ
P_SampleCov(1,1) = 2;                    % ��ʼ�����ֲ��ķ���

ParticleWeight    = zeros(1,NumParticle);     % ��ʼ��Ȩ��
x_SamplePart_temp = zeros(1,NumParticle);     % �м����
x_SampleParticle  = zeros(NumParticle,SimuTimes);

% Intinialization particle, prior distirbution p(x_0) 
for i = 1 : NumParticle
    % ��ʼ״̬���� x=0 ��ֵ������Ϊ sqrt(P) �ĸ�˹�ֲ�
    x_SampleParticle(i,1) = x_EstiState(1,1) + q + sqrt(P_SampleCov(1,1)) * randn; 
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
% xhatPartArr = [xhatPart]; %

f = @(x)3*sin(2*x) - exp(x);
h = @(x)x;

% ����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient
bino_fir = zeros(1,SimuTimes);       % ΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:NumParticle
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

%%
% diff_X_real ��ʾkʱ��״̬��΢��
diff_X_real = 0;

tic

%% ����ʵ��״̬ calculate real state and measurement
for k = 2 : SimuTimes
    diff_X_real = f(X_RealState(1,k-1)) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i) * X_RealState(1,k+1-i);
    end
    X_RealState(1,k) = diff_X_real - rema;
    % k ʱ����ʵֵ
    Y_RealMeas(1,k) = h(X_RealState(1,k)) + V_noise(1,k);  % k ʱ�̹۲�ֵ

%% ����N������ 
 for i = 1 : NumParticle
     % Draw particle: x^i_k ~ p(x_k | x^i_k-1) state transform function
     % ������� Num_particle ������
     x_SamplePart_temp(1,i) =  f(x_SampleParticle(i,k-1)) + q + sqrt(Q) * randn;
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + bino_fir(1,j)*x_SampleParticle(i,k+1-j);
         end
     x_SamplePart_temp(1,i) = x_SamplePart_temp(1,i) - temp;
     y_ParticleMeas = h(x_SamplePart_temp(1,i)) + r;     % ÿ�����Ӷ�Ӧ�Ĺ۲�ֵ
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
              x_SampleParticle(i,k) = x_SamplePart_temp(1,j);
              break;
          %else
          %    x_SampleParticle(i,k) = x_SampleParticle(i,k-1);
          end
      end
  end

%% ����ϵͳ״̬
  
x_EstiState(1,k) = mean(x_SampleParticle(:,k));

end


 FPF_SISO_TIME(1,kkk)    = toc;
 FPF_ERROR_norm1(1,kkk)  = norm((X_RealState - x_EstiState),1);
 FPF_ERROR_norm2(1,kkk)  = norm((X_RealState - x_EstiState),2);
 FPF_ERROR_RMSE(1,kkk) = sqrt(sum((X_RealState - x_EstiState).^2)/50);

 end

 FPF100_SISO_ERROR_TIME(1,:) = FPF_ERROR_norm1(1,:);
 FPF100_SISO_ERROR_TIME(2,:) = FPF_ERROR_norm2(1,:);
 FPF100_SISO_ERROR_TIME(3,:) = FPF_SISO_TIME(1,:);
 FPF100_ERROR_norm1_average  = sum(FPF100_SISO_ERROR_TIME(1,:))/50
 FPF100_ERROR_norm2_average  = sum(FPF100_SISO_ERROR_TIME(2,:))/50
 FPF100_SISO_TIME_average    = sum(FPF100_SISO_ERROR_TIME(3,:))/50
  save FPF_SISO_ERROE_TIME1 FPF100_SISO_ERROR_TIME FPF100_ERROR_norm1_average ...
       FPF100_ERROR_norm2_average FPF100_SISO_TIME_average




