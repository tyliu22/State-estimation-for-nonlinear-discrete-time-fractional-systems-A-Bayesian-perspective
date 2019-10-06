%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOSystemDataSet\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));
% SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));
Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%���沽��
N = SystemParameter{1,1}{1};     % ���沽��
q = SystemParameter{1,1}{2};     % ϵͳ������ֵ
r = SystemParameter{1,1}{3};     % ����������ֵ
Q = SystemParameter{1,1}{4};     % ϵͳ��������
R = SystemParameter{1,1}{5};     % ������������
alpha = SystemParameter{1,1}{6}; % ϵͳ�״�

%GL�����¶̼���ԭ��ĳ���
L = N+1;

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

% ϵͳ�������������
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
   
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------�����׿������˲������ܲ���------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_con = sqrt(3);

X_esti = zeros(3,N);          %״̬���Ź���ֵ
X_esti(:,1) = [0.1 0.1 0.1]'; %��ʼ״̬����
P_xesti = cell(1,N);         %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = [10,0,0;0,10,0;0,0,10];   %��ʼԤ�ⷽ����
P_xesti{1,1} = P_pred_0;             %��ʼ���Ʒ�����

for k=2:1:N
      %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = f(X_esti(:,k-1));
            %��������
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
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_xesti{1,k-1})';

            %��������
            rema_P = [0,0,0;0,0,0;0,0,0];
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_xesti{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + gamma{1,i}*P_xesti{1,k+1-i}*gamma{1,i}';
                end
            end

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun1 = f(X_esti(:,k-1)+h_con*S_chol(:,1)) - f(X_esti(:,k-1)-h_con*S_chol(:,1));
        temp_fun2 = f(X_esti(:,k-1)+h_con*S_chol(:,2)) - f(X_esti(:,k-1)-h_con*S_chol(:,2));
        temp_fun3 = f(X_esti(:,k-1)+h_con*S_chol(:,3)) - f(X_esti(:,k-1)-h_con*S_chol(:,3));
        
        temp_fun = temp_fun1 * temp_fun1' + temp_fun2 * temp_fun2' + temp_fun3 * temp_fun3';
        
        temp_initial = ( 0.5*temp_fun1*S_chol(:,1)'+ ...
                         0.5*temp_fun2*S_chol(:,2)'+ ...
                         0.5*temp_fun3*S_chol(:,3)' ...
                         )*gamma{1,2}';
        
        temp = 1/(4*h_con^2) * temp_fun + rema_P + ...
                  - 1/h_con * temp_initial - 1/h_con * temp_initial';
        P_pred = temp + Q;
        
        %��״̬Ԥ�����Э����������cholsky�ֽ�
        S_chol_pred = chol(P_pred)';
        
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        
        temp_fun1 = h(X_pre+h_con*S_chol_pred(:,1)) - h(X_pre-h_con*S_chol_pred(:,1));
        temp_fun2 = h(X_pre+h_con*S_chol_pred(:,2)) - h(X_pre-h_con*S_chol_pred(:,2));
        temp_fun3 = h(X_pre+h_con*S_chol_pred(:,3)) - h(X_pre-h_con*S_chol_pred(:,3));
        
        P_zpred = 1/(4*h_con^2) * ( temp_fun1 * temp_fun1' + temp_fun2 * ...
                  temp_fun2' + temp_fun3 * temp_fun3' ) + R;

        %Ԥ����Э�������:P_xzpred
        P_xzpred = 1/(2*h_con) * ( S_chol_pred(:,1)*temp_fun1' ...
                  + S_chol_pred(:,2)*temp_fun2' + S_chol_pred(:,3)*temp_fun3');

        %���㿨��������:Kk(2*1)
        Kk = P_xzpred/P_zpred;

        %״̬����
        X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_xesti{1,k} = P_pred - Kk*P_zpred*Kk';
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

FCDKF_X_esti = X_esti;

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOEstimatedState\';
save(strcat(path,'FCDKF_EstimatedState','.mat'),'FCDKF_X_esti');





