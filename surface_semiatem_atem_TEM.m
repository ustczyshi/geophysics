%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平电偶源激励，地面或空中观测
%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;

%   for k2 = 1:5
% d = [0 25 30 35 50];
%%
mu_0 = 4*pi*1e-7;
rou = 1;
%% 地面
r=50; % 收发水平偏移距
a = 10; % 圆形回线半径
% 半航空收发高度参数
z =0;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度


%% 全航空收发高度参数
%% 仿真场景选择
%{
disp('************************************************');
disp('仿真模型选择：');
disp('0：均匀大地地面TEM----------');
disp('1：均匀大地全航空TEM----------');
disp('2：均匀大地半航空TEM仿真----------');
tig = input('请选择：');
switch tig
    case 0
        r  = input('收发水平距：');
%         a = input('圆形回线半径 ：');
        h = 0;
        z = 0;
    case 1  %  均匀大地全航空TEM
%         a = input('圆形回线半径 ：');
        h = input('发射源距地面高度：');
        z = input('观测线圈距地面高度：'); 
    case 2   %  均匀大地半航空TEM
        h = 0;
        z = input('观测线圈距地面高度：'); 
    case -1
        h = 0;
        z = 0;
        r = 100;
    otherwise 
         error('输入选择错误！默认仿真地面TEM ');
         tig = input('请选择：');
        return;
end
disp('************************************************');
%}
%% 源激励模型选择: 1--垂直磁偶极源；0：圆形回线源

% fig = 1;% 垂直磁偶极源
 fig = 0;% 圆形回线
%% 
G_S=load ('G_S.txt')';
m2 = 1:length(G_S);
fs = 1e7;
t = 1/fs:1/fs:1e-2;% 时间区间
% t = logspace(-8,1,10000);
%%  z轴变量
h_z1_t=zeros(1,length(t));% 负脉冲响应时域
h_z2_t=h_z1_t;%负阶跃
h01_z1_t=zeros(1,length(t));% 负脉冲响应时域
h01_z2_t=h_z1_t;%负阶跃
%%  x轴变量
h_x1_t=zeros(1,length(t));
h_x2_t=h_x1_t;
h01_x1_t=zeros(1,length(t));% 负脉冲响应时域
h01_x2_t=h_x1_t;%负阶跃
%%  y方向
h_y1_t=zeros(1,length(t));
h_y2_t=h_y1_t;
h01_y1_t=zeros(1,length(t));% 负脉冲响应时域
h01_y2_t=h_y1_t;%负阶跃
%--------------------------------------------------------------------------
%第一步：读取已经存储的滤波器系数,表示为行向量；
%%  y方向计算
load J0_Gupt.txt;       
J_zero = J0_Gupt( :, 3)'; % 快速汉克尔变换滤波系数
delta = J0_Gupt( :, 2)'; %  采样点的横坐标偏移量
%% x-z方向计算
 load J1_Gupt.txt;       
J_1 = J1_Gupt( :, 3)'; % 快速汉克尔变换滤波系数
delta_1 = J1_Gupt( :, 2)'; %  采样点的横坐标偏移量
%--------------------------------------------------------------------------
%计算lambda，并将lambda和frequency扩展成二维矩阵-----针对垂直方向汉克尔变换
if fig == 1; % 水平磁偶极源
    lambda = (1./r) .*exp(delta); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda，针对J0
    lambda_r = (1./r) .*exp(delta_1); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda,针对J1汉克尔变换
end
%{
if fig ==0 ;% 圆形回线，垂直方向磁场汉克尔变换的lambda值 
    lambda = (1./a) .*exp(delta_1); 
end
%}
%--------------------------------------------------------------------------

    for ii=1:length(t)
        freq = (log(2)*1i/(t(ii)*2*pi))*m2;
        %--------------------------------------------------------------------------
        %根据递推公式计算空气-地面的反射系数
        if fig == 1;
            [lambda_Array,frequency_Array] = Array_trans(lambda,freq); % 针对垂直方向汉克尔变换
            [lambdar_Array,frequencyr_Array] = Array_trans(lambda_r,freq); % 针对水平方向汉克尔变换   
            r_TE=calculate_r_TE(lambda,freq); % 针对垂直方向汉克尔变换
            r_TEr=calculate_r_TE(lambda_r,freq); % 针对垂直方向汉克尔变换
        end
        if fig ==0;
            [lambda_Array,frequency_Array] = Array_trans(lambda,freq); % 针对垂直方向汉克尔变换
            r_TE=calculate_r_TE(lambda,freq); % 针对垂直方向汉克尔变换
        end
        %% --------------------------------------------------------------------------
        %{
        %计算磁场的垂直分量，并用快速汉克尔变换求积分的数值解。  
        %垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置
        % 垂直磁偶极源z轴磁场的核函数 u_0 = lambda,z=h=0,地面发射地面接收;
    %     sum = (1+r_TE).*lambda_Array.^2;  % 正阶跃响应频域核函数
    %     sum_zeros = (1).*lambda_Array.^2; %% 对应直流成分，此时的核函数中的r_TE=0
    %     H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%正阶跃响应频域
    %     H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%正阶跃响应零频响应
    %     h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
    %     h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%负阶跃响应时域
        %}
        %% ---------------------------------------------------------------------------
        if fig == 1; % 垂直磁偶极源
            %% -----------------------------------半航空tem 垂直方向----------------------------------------
            % 全航空atem垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置.对于全航空,z=-50,h=-100米。
            sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.^2;
            sum_zeros = (1).*lambda_Array.^2.*exp(-lambda_Array*(z+h)); %% 对应直流成分，此时的核函数中的r_TE=0
            H_vertical = 1./(4*pi) *  Fast_Hankel(r,sum,J_zero);%正阶跃响应频域
            H_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sum_zeros,J_zero);%正阶跃响应零频响应
            h_z1_t(ii) = -GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_z2_t(ii) = +GS_Trans(t(ii),H_vertical_zeros,freq,G_S)-GS_Trans(t(ii),H_vertical,freq,G_S);%负阶跃响应时域
            h01_z1_t(ii) = GS_Trans2(t(ii),H_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_z2_t(ii) = GS_Trans(t(ii),H_vertical,freq,G_S);%正阶跃响应时域
            %% -----------------------------------半航空tem 径向----------------------------------------
            sumr = (exp(-lambdar_Array*(z+h))-r_TEr.*exp(lambdar_Array*(z-h))).*lambdar_Array.^2;
            sumr_zeros = (1).*lambdar_Array.^2.*exp(-lambdar_Array*(z+h)); %% 对应直流成分，此时的核函数中的r_TE=0
            Hr_vertical = 1./(4*pi) *  Fast_Hankel(r,sumr,J_1);%正阶跃响应频域
            Hr_vertical_zeros = 1./(4*pi) *  Fast_Hankel(r,sumr_zeros,J_1);%正阶跃响应零频响应
            h_r1_t(ii) = -GS_Trans2(t(ii),Hr_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_r2_t(ii) = +GS_Trans(t(ii),Hr_vertical_zeros,freq,G_S)-GS_Trans(t(ii),Hr_vertical,freq,G_S);%负阶跃响应时域
            h01_r1_t(ii) = GS_Trans2(t(ii),Hr_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_r2_t(ii) = GS_Trans(t(ii),Hr_vertical,freq,G_S);%正阶跃响应时域
            %% -----------------------------------半航空tem phi方向电场-----------------------------------
            sump = (exp(-lambdar_Array*(z+h))+r_TEr.*exp(lambdar_Array*(z-h))).*lambdar_Array;
            sump_zeros = (1).*lambdar_Array.*exp(-lambdar_Array*(z+h)); %% 对应直流成分，此时的核函数中的r_TE=0
            Hp_vertical = -1i*2*pi*mu_0*freq'./(4*pi) .* Fast_Hankel(r,sump,J_1);%正阶跃响应频域
            Hp_vertical_zeros = -1i*2*pi*mu_0*freq'./(4*pi) .* Fast_Hankel(r,sump_zeros,J_1);%正阶跃响应零频响应
            h_p1_t(ii) = -GS_Trans2(t(ii),Hp_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_p2_t(ii) = +GS_Trans(t(ii),Hp_vertical_zeros,freq,G_S)-GS_Trans(t(ii),Hp_vertical,freq,G_S);%负阶跃响应时域
            h01_p1_t(ii) = GS_Trans2(t(ii),Hp_vertical,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_p2_t(ii) = GS_Trans(t(ii),Hp_vertical,freq,G_S);%正阶跃响应时域
        end
        if fig == 0; % 圆形回线
            %% ***************************垂直方向 z 轴*********************************
            sum = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.*besselj(0,lambda_Array.*r);
            sum_zeros = exp(-lambda_Array*(z+h)).*lambda_Array.*besselj(0,lambda_Array.*r);%% 对应直流成分，此时的核函数中的r_TE=0
            H_z = 0.5*a * Fast_Hankel(a,sum,J_1);%正阶跃响应频域
            H_z_zeros = 0.5*a *  Fast_Hankel(a,sum_zeros,J_1);
            h_z1_t(ii) = -GS_Trans2(t(ii),H_z,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_z2_t(ii) = +GS_Trans(t(ii),H_z_zeros,freq,G_S)-GS_Trans(t(ii),H_z,freq,G_S);%负阶跃响应时域
            h01_z1_t(ii) = GS_Trans2(t(ii),H_z,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_z2_t(ii) = GS_Trans(t(ii),H_z,freq,G_S);%正阶跃响应时域
            %% **************************径向 r 磁场**********************************************
            sumr = (exp(-lambda_Array*(z+h))-r_TE.*exp(lambda_Array*(z-h))).*lambda_Array.*besselj(1,lambda_Array.*r);
            sumr_zeros = exp(-lambda_Array*(z+h)).*lambda_Array.*besselj(1,lambda_Array.*r);%% 对应直流成分，此时的核函数中的r_TE=0
            H_r = 0.5*a * Fast_Hankel(a,sumr,J_1);%正阶跃响应频域
            H_r_zeros = 0.5*a *  Fast_Hankel(a,sumr_zeros,J_1);
            h_r1_t(ii) = -GS_Trans2(t(ii),H_r,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_r2_t(ii) = +GS_Trans(t(ii),H_r_zeros,freq,G_S)-GS_Trans(t(ii),H_r,freq,G_S);%负阶跃响应时域
            h01_r1_t(ii) = GS_Trans2(t(ii),H_r,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_r2_t(ii) = GS_Trans(t(ii),H_r,freq,G_S);%正阶跃响应时域
%             %% **************************phi方向电场**********************************************
            sump = (exp(-lambda_Array*(z+h))+r_TE.*exp(lambda_Array*(z-h))).*besselj(1,lambda_Array.*r);
            sump_zeros = exp(-lambda_Array*(z+h)).*besselj(1,lambda_Array.*r);%% 对应直流成分，此时的核函数中的r_TE=0
            H_p = -1i*2*pi*mu_0*freq'.*0.5*a .* Fast_Hankel(a,sump,J_1);%正阶跃响应频域
            H_p_zeros = -1i*2*pi*mu_0*freq'.*0.5*a .*  Fast_Hankel(a,sump_zeros,J_1);
            h_p1_t(ii) = -GS_Trans2(t(ii),H_p,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %负脉冲响应时域
            h_p2_t(ii) = +GS_Trans(t(ii),H_p_zeros,freq,G_S)-GS_Trans(t(ii),H_p,freq,G_S);%负阶跃响应时域
            h01_p1_t(ii) = GS_Trans2(t(ii),H_p,G_S);%GS_Trans2(t(ii),H_vertical_zeros); %正脉冲响应时域
            h01_p2_t(ii) = GS_Trans(t(ii),H_p,freq,G_S);%正阶跃响应时域
        end
    end
 %% 保存正脉冲响应和正阶跃响应
hzimpulse01 = h01_z1_t;
hzstep01 = h01_z2_t;
hrimpulse01 = h01_r1_t;
hrstep01 = h01_r2_t;
epimpulse01 = h01_p1_t;
epstep01 = h01_p2_t;
 %% 保存负脉冲响应和负阶跃响应
hzimpulse10 = h_z1_t;
hzstep10 = h_z2_t;
hrimpulse10 = h_r1_t;
hrstep10 = h_r2_t;
epimpulse10 = h_p1_t;
epstep10 = h_p2_t;
%% 保存生成数据
if fig ==1 
    save(['atem_magnetic_dipole_response_h' num2str(h) 'z' num2str(z) '_r' num2str(r) '_rou' num2str(rou) '.mat'],...
        'hzimpulse10','hzstep10','hrimpulse10','hrstep10','epimpulse10','epstep10',...
        'hzimpulse01','hzstep01','hrimpulse01','hrstep01','epimpulse01','epstep01');
end
if fig ==0 
    %% h0z0:地面发射地面接收--中心回线
    save(['atem_circle_loop_response_h' num2str(h) 'z' num2str(z) '_r' num2str(r) '_rou' num2str(rou) '.mat'],...
        'hzimpulse10','hzstep10','hrimpulse10','hrstep10','epimpulse10','epstep10',...
        'hzimpulse01','hzstep01','hrimpulse01','hrstep01','epimpulse01','epstep01');
end
%% --------------------------------------------------------------------------
% 解析解
I = 1;
rou = 100;
 [step_H01,step_H10,impluse_H10] = Central_Circle_Loop_Step_Response(mu_0,rou,I,a,t);
 
%% 负脉冲和负阶跃
%作图。
figure;
loglog(t.*10^3,abs(h_z1_t),'r','Linewidth',1);
hold on
loglog(t.*10^3,(h_z2_t),'b','Linewidth',1);
hold on 
loglog(t.*1000,abs(diff([0 step_H10])*fs),'r--','Linewidth',1);
hold on;
loglog(t.*1000,step_H10,'b--','Linewidth',1);
hold on;
grid on;
legend('数值解正脉冲响应','数值解正阶跃响应','解析解正脉冲响应','解析解正阶跃响应');
title('负阶跃圆形回线的垂直磁场')
xlabel('时间（ms）')
ylabel('Hz/(A/m)');
%% 正阶跃
%作图。
figure;
loglog(t.*10^3,abs(h01_z1_t),'r','Linewidth',1);
hold on
loglog(t.*10^3,(h01_z2_t),'b','Linewidth',1);
hold on 
loglog(t(1:end-1).*1000,abs(diff(step_H01)*fs),'r--','Linewidth',1);
hold on;
loglog(t.*1000,step_H01,'b--','Linewidth',1);
hold on;
grid on;
legend('数值解正脉冲响应','数值解正阶跃响应','解析解正脉冲响应','解析解正阶跃响应');
title('正阶跃圆形回线的垂直磁场');
xlabel('时间（ms）');
ylabel('Hz/(A/m)');
%% --------------------------------------------------------------------------
%作图。
figure;
loglog(t.*10^3,abs(h_r1_t),'r','Linewidth',2);
hold on;
loglog(t.*10^3,0.5*(abs(h_r1_t)-h_r1_t),'r--','Linewidth',2);
hold on;
loglog(t.*10^3,abs(h_r2_t),'b','Linewidth',2);
hold on ;
loglog(t.*10^3,0.5*(abs(h_r2_t)-h_r2_t),'b--','Linewidth',2);
hold on ;
grid on;
title('瞬断圆形回线的径向磁场');
xlabel('时间（ms）');
ylabel('Hr/(A/m)');
%% --------------------------------------------------------------------------
%作图。
figure;
loglog(t.*10^3,abs(h_p1_t),'r','Linewidth',2);
hold on;
loglog(t.*10^3,0.5*(abs(h_p1_t)-h_p1_t),'r--','Linewidth',2);
hold on;
loglog(t.*10^3,abs(h_p2_t),'b','Linewidth',2);
hold on ;
loglog(t.*10^3,0.5*(abs(h_p2_t)-h_p2_t),'b--','Linewidth',2);
hold on ;
grid on;
title('瞬断圆形回线的phi电场');
xlabel('时间（ms）');
ylabel('Ep/(A/m)');
close all;clear all;
%  end

