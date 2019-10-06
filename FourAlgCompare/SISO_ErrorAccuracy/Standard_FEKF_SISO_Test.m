%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order EKF
%   Ŀ�ģ�FCDKF��FEKF�����ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save��load��ʹ�÷�����
% save(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'),'data');
% load(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'));

% SystemParameter = zeros(1,6);
% SystemParameter(1) = N;     % ���沽��
% SystemParameter(2) = q;     % ϵͳ������ֵ
% SystemParameter(3) = r;     % ����������ֵ
% SystemParameter(4) = Q;     % ϵͳ��������
% SystemParameter(5) = R;     % ������������
% SystemParameter(6) = alpha; % ϵͳ�״�

clc
clear

% ָ����ǰ·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';
SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));

%���沽��
N = SystemParameter(1);     % ���沽��
q = SystemParameter(2);     % ϵͳ������ֵ
r = SystemParameter(3);     % ����������ֵ
Q = SystemParameter(4);     % ϵͳ��������
R = SystemParameter(5);     % ������������
alpha = SystemParameter(6); % ϵͳ�״�

Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%GL�����¶̼���ԭ��ĳ���
L = N+1;

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ alpha ʱGL�����µ�ϵ��
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%ϵͳ��������
I = eye(1,1);                %���ɵ�λ��

% ϵͳ�������������
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------��������չ�������˲������ܲ���---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti = zeros(1,N);        %״̬���Ź���ֵ
P_xesti = zeros(1,N);       %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = 20;              %��ʼԤ�ⷽ����
P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����

for k=2:1:N
  %�������˲�
      %״̬Ԥ��:X_pre
        diff_X_esti = f(X_esti(1,k-1));
            %��������
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
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
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
        
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;
        
        %���㿨��������:Kk(2*1)
        H = 1;
        Kk = P_xpred*H'/(H*P_xpred*H' + R);
        
        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
        
        %���Ʒ���������
        P_xesti(1,k) = (I-Kk*H)*P_xpred;
end

% %������������ͼ
% k = 1:1:N;
% 
% LineWidth = 1.5;
% 
% figure;
% plot(k,X_real(1,:),'r',k,X_esti(1,:),'b--','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 N -6 6 ])
% ylabel('x','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% % ����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');

FEKF_X_esti = X_esti;

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
save(strcat(path,'FEKF_EstimatedState','.mat'),'FEKF_X_esti');
