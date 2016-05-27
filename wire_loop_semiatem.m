%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平电偶源激励，地面或空中观测

%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%%
u0 = 4*pi*1e-7;
load parameters.txt;
sigma1 = parameters(1,2);%第一层的电导率
rou = 1./sigma1;
%% 发射机参数地面
x = 100; % 收发水平偏移距，沿y轴
y = 100;
L= 1; % 发射线缆长度，沿x轴
I = 1; % 发射电流
%%  半航空收发高度参数
z =0;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度

%% 采样率和观测时间段设置
fs = 1e5;% 采样率
dt = 1./fs;
t = 1/fs:1/fs:4e-2;% 时间区间
%% 电场分量
%  x轴变量
% E_x_impulse=zeros(1,length(t));% 负脉冲响应时域
% E_x_step = E_x_impulse;%负阶跃
% E_x_impulse1=zeros(1,length(t));% 正脉冲响应时域
% E_x_step1 = E_x_impulse1;%正阶跃
% %  y方向
% E_y_impulse=zeros(1,length(t));% 负脉冲响应时域
% E_y_step = E_y_impulse;%负阶跃
% E_y_impulse1=zeros(1,length(t));% 正脉冲响应时域
% E_y_step1 =E_y_impulse1;%正阶跃
% %% 磁场分量
% %  z轴变量
% h_z_impulse=zeros(1,length(t));% 负脉冲响应时域
% h_z_step = h_z_impulse;%负阶跃
% h_z_impulse1=zeros(1,length(t));% 正脉冲响应时域
% h_z_step1 = h_z_impulse1;%正阶跃
% %  x轴变量
% h_x_impulse=zeros(1,length(t));% 负脉冲响应时域
% h_x_step = h_x_impulse;%负阶跃
% h_x_impulse1=zeros(1,length(t));% 正脉冲响应时域
% h_x_step1 = h_x_impulse1;%正阶跃
% %  y方向
% h_y_impulse=zeros(1,length(t));% 负脉冲响应时域
% h_y_step = h_y_impulse;%负阶跃
% h_y_impulse1=zeros(1,length(t));% 正脉冲响应时域
% h_y_step1 = h_y_impulse1;%正阶跃
%%
% [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t)
[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
% save('horizontal_electrical_dipole_impulse_shuzhijie','hx_1_impulse','hy_1_impulse','hz_1_impulse');
%% z轴阶跃
%{
%    作图。
hz_10 = hz_01(end) - hz_01;
figure;
plot(t.*10^3,u0.*hz_01,'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*hz_10,'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x轴阶跃
%作图。
% hx_10 = hx_01(end) - hx_01;
figure;
plot(t.*10^3,u0.*(hx_01),'r','Linewidth',1);
hold on
plot(t.*10^3,-u0.*(hx_01),'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y 轴阶跃
hy_10 = hy_01(end) - hy_01;
figure;
plot(t.*10^3,u0.*(hy_01),'r','Linewidth',1);
hold on
plot(t.*10^3,u0.*(hy_10),'b','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%}
%% z 轴脉冲
load horizontal_electrical_dipole_impulse_jiexijie.mat;
figure;
loglog(t(1:end).*10^3,abs(real(ex_01)),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ex01),'k:','Linewidth',2);
grid on;
legend('数值解ex\_01','解析解step\_ex01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
loglog(t(1:end).*10^3,abs(real(ey_01)),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
grid on;
legend('数值解ey\_01','解析解step\_ey01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
% 误差分析
% ex_error = abs(abs(ex_01)-abs(step_ex01))./abs(step_ex01);
% ey_error = abs(abs(ey_01)-abs(step_ey01))./abs(step_ey01);
% figure;
% loglog(t.*1e3,ex_error.*100,'r','linewidth',2);
% hold on;
% loglog(t.*1e3,ey_error.*100,'k:','linewidth',2);
% grid on;
% legend('ex','ey');
% title('数值解和解析解的误差分析');
% xlabel('Time/(ms)')
% ylabel('error/(%)');
% 脉冲响应
% figure;
% plot(t(1:end).*10^3,ex_impulse,'r','Linewidth',2);
% hold on
% plot(t(1:end-1).*10^3,diff(step_ex01)./dt,'k:','Linewidth',2);
% grid on;
% legend('数值解ex\_01','解析解step\_ex01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
% xlabel('Time/(ms)')
% ylabel('Ex/(V/m)');
% figure;
% plot(t(1:end).*10^3,abs(ey_01),'r','Linewidth',2);
% hold on
% plot(t(1:end).*10^3,abs(step_ey01),'k:','Linewidth',2);
% grid on;
% legend('数值解ey\_01','解析解step\_ey01');
% title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
% xlabel('Time/(ms)')
% ylabel('Ey/(V/m)');
%
%%
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hz),'k:','Linewidth',2);
grid on;
legend('数值解hz\_1\_impulse','数值解diff(hz\_01)./dt','解析解impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hx_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hx_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hx),'k:','Linewidth',2);
grid on;
legend('数值解hx\_1\_impulse','数值解diff(hx\_01)./dt','解析解impulse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
figure;
loglog(t(1:end).*10^3,u0.*abs(hy_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hy_01)./dt),'b','Linewidth',2);
grid on;
loglog(t(1:end).*10^3,u0.*abs(impulse_hy),'k:','Linewidth',2);
grid on;
legend('数值解hy\_1\_impulse','数值解diff(hy\_01)./dt','解析解impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
% 误差分析
hz_error = abs(abs(hz_1_impulse)-abs(impulse_hz))./abs(impulse_hz);
hx_error = abs(abs(hx_1_impulse)-abs(impulse_hx))./abs(impulse_hx);
hy_error = abs(abs(hy_1_impulse)-abs(impulse_hy))./abs(impulse_hy);
figure;
loglog(t.*1e3,hz_error.*100,'r','linewidth',2);
hold on;
loglog(t.*1e3,hx_error.*100,'b','linewidth',2);
hold on;
loglog(t.*1e3,hy_error.*100,'k:','linewidth',2);
grid on;
legend('hz','hx','hy');
title('数值解和解析解的误差分析');
xlabel('Time/(ms)')
ylabel('error/(%)');

%% save data

