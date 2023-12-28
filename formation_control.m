clear
clc
% ����ʱ��
e = 0.02;
% ������
d = 0.1;
% ���Ӱ뾶
rw = 0.1;
% ����뾶
Rw = 0.5; 
% ����
m = 1;
% ת������
I = 0.1;
% δ֪�Ŷ�
td = [0.0; 0.0; 0.00];
% ����
R = [1 0
     0 1];
% ���س���
km = [0.5 0
      0   0.5];
% ���������ԽǾ���
Jm = [0.002 0
      0     0.002];
% reduction gear ���ٳ���
r = [0.5 0
     0   0.5];
% ���
La = [0.00021 0
      0       0.00021];
% ���綯��
kb = [0.5 0
      0   0.5];
% ˥��
Bm = [0.001 0
      0     0.001];%�ԽǾ���

% Ħ����
f = [0; 0; 0];

% ����ʱ���ܳ�
N = 6000;
% ref�ĳ�ʼ״̬ q = [x y theta]��������һ��Բ��
qr(:,1) = [cos(0/100*2*pi); sin(0/100*2*pi);0];
vqr(:,1) = [0; 0];
for i=1:N
    S = [cos(qr(3,i)) -d*sin(qr(3,i))
         sin(qr(3,i)) d*cos(qr(3,i))
         0 1];
    vqr(:,i) = [2*pi/100; 2*pi/100];
%     vqr(:,i) = [0.01; 0];
    qr(:,i+1) = qr(:,i) + S * vqr(:,i) * e;
end
plot(qr(1,:),qr(2,:))
axis equal
hold on 
% leader�ĳ�ʼ״̬
qi(:,1) = [0; 0; 0];
vqi(:,1) = [0; 0];
eic(:,1) = [0; 0];
% �ο����ٶ�
vic(:,1) = [0; 0];
ui(:,1) = [0; 0];

% ���Ʋ���
ki1= 1; ki2= 1; ki3= 1;
for k=1:N
    % leader��ǰ�Ƕ�
    theta = qi(3,k);
    % leader��ǰ���ٶ�
    theta_dot = vqi(2,k);
    
    % �ſɱȾ��������ֽ��ٶ���ѿ�������ϵ�Ĺ�ϵ    
    Jq     = (rw/2) * [cos(theta) cos(theta)
                       sin(theta) sin(theta)
                       1/Rw       -1/Rw];
    Jq_dot = (rw/2) * [-sin(theta)*theta_dot -sin(theta)*theta_dot
                        cos(theta)*theta_dot  cos(theta)*theta_dot
                        0           0];
    
    S = [cos(theta) -d*sin(theta)
         sin(theta)  d*cos(theta)
         0           1];
    % ���Ծ���     
    M = [m               0              m*d*sin(theta)
         0               m             -m*d*cos(theta)
         m*d*sin(theta) -m*d*cos(theta) I];
     
    % �������������������
    V = [m * d * theta_dot^(2) * cos(theta)
         m * d * theta_dot^(2) * sin(theta)
         0]';
    % ת������ϵ��
    B = (1 / rw) *[cos(theta) cos(theta)
                   sin(theta) sin(theta)
                   Rw        -Rw];
    % תΪ���ٶȣ����ٶ������صķ��̣�M*v_dot + Vm*v+F_+td_
    M_ = S'*M*S;
    S_dot = [-sin(theta)*theta_dot -d*cos(theta)*theta_dot
              cos(theta)*theta_dot -d*sin(theta)*theta_dot
              0           0]; 
    Vm_ = S' * (M * S_dot + V * S);
    B_ = S' * B;
    td_ = S' * td;
    f_ = S' * f;
    %��ϵ�������
    Jq_plus = (1 / rw) * [cos(theta) sin(theta)  Rw
                          cos(theta) sin(theta) -Rw];
    Jq_plusdot = (theta_dot / rw) * [-sin(theta) cos(theta) 0
                                     -sin(theta) cos(theta) 0];
                                 
    Dq = R * km^-1 * (Jm * r^-1 * Jq_plus * S + r * B_^-1 * M_);
    Nq = R * km^-1 * (Jm * r^-1 * Jq_plus * S_dot + Jm * r^-1 * Jq_plusdot * S ...
         + Bm * r^-1 * Jq_plus * S + r * B_^-1 * Vm_) + kb * r^-1 * Jq_plus * S;
    Hq = R * km^-1*r*B_^-1*(f_+ td_);
    % vl��δ��ģ��̬��vL=La*Ia_dot
    vL=[0;0];
    % ����ƫ��ֵ 
    Tei = [ cos(qi(3,k)) sin(qi(3,k)) 0
           -sin(qi(3,k)) cos(qi(3,k)) 0
            0            0            1];
    ei = Tei * (qr(:,k)-qi(:,k));
    
    vr = vqr(:,k);
    vic(:,k+1) = [vqr(1,k) * cos(ei(3)) + ki1 * ei(1)
                  vqr(2,k) + ki2 * vqr(1,k) * ei(2) + ki3 * vqr(1,k) * sin(ei(3))];
 
    eic(:,k) = vic(:,k) - vqi(:,k);  
    Ki4 = 25;
    ui(:,k) = Ki4 * eic(:,k); 
    vqi(:,k+1) = vqi(:,k) + ui(:,k) * e;
    
    qi(:,k+1) = qi(:,k) + S*vqi(:,k+1)*e;
%     qi(:,k+1) = qi(:,k) + S*vic(:,k+1)*e;
    plot(qi(1,:),qi(2,:)) 
    pause(0.01)
end 


plot(qi(1,:),qi(2,:))




