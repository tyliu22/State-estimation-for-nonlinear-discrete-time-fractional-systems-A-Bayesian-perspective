%*************************************************************************%
%   �����������˲����渴��
%   ���ģ�     fractional order PF
%   Ŀ�ģ������������˲��㷨����
%         �ԷǸ�˹ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й���
%
%   ��ע�����ȶ��ֲ�
%          
%*************************************************************************%

clc;
clear all;

LineWidth = 1.5;

SimuTimes = 50;    % ����ʱ��
NumParticle = 100;  % ���Ӹ���

%ϵͳ��������
I = eye(1,1);               % ���ɵ�λ��

%����
P_alpha_W = 1.6; % alpha�ȶ��ֲ� ����ָ��
P_beta_W  = 0.8; % alpha�ȶ��ֲ� �߶Ȳ���
P_gamma_W = 0.1; % alpha�ȶ��ֲ� ƫб����
P_delta_W = 0;   % alpha�ȶ��ֲ� λ�ò���

P_alpha_V = 1.8; % alpha�ȶ��ֲ� ����ָ��
P_beta_V  = 0.9; % alpha�ȶ��ֲ� �߶Ȳ���
P_gamma_V = 0.2; % alpha�ȶ��ֲ� ƫб����
P_delta_V = 0;   % alpha�ȶ��ֲ� λ�ò���

W_noise = (stblrnd(P_alpha_W,P_beta_W,P_gamma_W,P_delta_W,SimuTimes,1))';  % ϵͳ����         % ϵͳ����
V_noise = (stblrnd(P_alpha_V,P_beta_V,P_gamma_V,P_delta_V,SimuTimes,1))';  % ��������

X_RealState = zeros(1,SimuTimes); % ϵͳ״̬��ʵֵ ��ʼֵ0
Y_RealMeas = zeros(1,SimuTimes);  % ϵͳ״̬��ʵֵ ��ʼֵ0
Y_RealMeas(1,1) = X_RealState(1,1) + ...
                  (stblrnd(P_alpha_V,P_beta_V,P_gamma_V,P_delta_V,1,1))';

% P_SampleCov = zeros(1,SimuTimes);        % ��������
x_EstiState = zeros(1,SimuTimes);          % ״̬����ֵ
% P_EstiState = zeros(1,SimuTimes);        % ״̬����
% P_SampleCov(1,1) = 2;                    % ��ʼ�����ֲ��ķ���

ParticleWeight    = zeros(SimuTimes,NumParticle);     % ��ʼ��Ȩ��
x_SamplePart_temp = zeros(SimuTimes,NumParticle);     % �м����
x_SampleParticle  = zeros(NumParticle,SimuTimes);

% Intinialization particle, prior distirbution p(x_0) 
for i = 1 : NumParticle
    % ��ʼ״̬���� x=0 ��ֵ������Ϊ sqrt(P) �ĸ�˹�ֲ�
    x_SampleParticle(i,1) = x_EstiState(1,1) + ...
                            stblrnd(P_alpha_W,P_beta_W,P_gamma_W,P_delta_W,1,1); 
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
     x_SamplePart_temp(k,i) = f(x_SampleParticle(i,k-1)) + ...
                              (stblrnd(P_alpha_W,P_beta_W,P_gamma_W,P_delta_W,1,1))';
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + bino_fir(1,j)*x_SampleParticle(i,k+1-j);
         end
     x_SamplePart_temp(k,i) = x_SamplePart_temp(k,i) - temp;
     y_ParticleMeas = h(x_SamplePart_temp(k,i));     % ÿ�����Ӷ�Ӧ�Ĺ۲�ֵ
     % ErrorMeas = Y_RealMeas(1,k) - y_ParticleMeas;   % ����ʵ�۲�֮�����Ȼ
     % Draw weight: w^i_k ~ p(z_k | x^i_k) measurement transform function
     % ����Ȩֵ������������й�
     %ParticleWeight(1,i) = h(x_SamplePart_temp(1,i)) + r + sqrt(R) * randn;
     ParticleWeight(k,i) = stblpdf( Y_RealMeas(1,k),P_alpha_V,P_beta_V,P_gamma_V,y_ParticleMeas,'quick');
     % ÿ�����ӵ���Ȼ�����ƶ�
 end

%%
 % Ȩֵ��һ��
weight_sum = sum(ParticleWeight(k,:));
for i = 1 : NumParticle
    ParticleWeight(k,i) = ParticleWeight(k,i) / weight_sum;  % ��һ�����Ȩֵ q
end

 % ����Ȩֵ���²��� ����ز����㷨
 qtempsum = zeros(1,NumParticle);
  qtempsum(1,1) = ParticleWeight(k,1);
 for i = 2 : 1 : NumParticle
    qtempsum(1,i) = qtempsum(1,i-1) + ParticleWeight(k,i);
 end
 
  for i = 1 : NumParticle
      UniRandom = rand; % �������ȷֲ������
      for j = 1 : NumParticle
          % �ۼ�Ȩֵ
          %qtempsum = qtempsum + ParticleWeight(1,j);
          if qtempsum(1,j) >= UniRandom
              x_SampleParticle(i,k) = x_SamplePart_temp(k,j);
              break;
          %else
          %    x_SampleParticle(i,k) = x_SampleParticle(i,k-1);
          end
      end
  end
 
  
%% ����ϵͳ״̬��Э����
x_EstiState(1,k) = mean(x_SampleParticle(:,k));
% P_EstiState(1,k) = sum((x_SampleParticle(:,k)-x_EstiState(1,k)).^2)/NumParticle;

end


%% estimation accuracy
t = 1 : SimuTimes;
figure;
plot(t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
hold on

for i = 1:NumParticle-1
    plot(t, x_SampleParticle(i,:),'c:');
    hold on
end
plot(t, x_SampleParticle(NumParticle,:),'c');
hold on
plot(t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
legend('real state','estimated state','particles');
%L2 = plot(t, x_SampleParticle(100,k),'y',t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
%Esitimated_state = legend('real state','estimated state','Location','best');
%set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
axis([0 50 -10 4]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('iteration times','FontSize',7); 
ylabel('state estimation','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)


%% estimation accuracy
t = 1 : SimuTimes;
figure;
plot(t, X_RealState, 'b', t, x_EstiState, 'r--','linewidth',LineWidth);
Esitimated_state = legend('real state','estimated state','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
axis([0 50 -10 4]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('iteration times','FontSize',7); 
ylabel('state estimation','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)




%% ���������
save X_RealState X_RealState
save x_EstiState x_EstiState
save x_SampleParticle x_SampleParticle




