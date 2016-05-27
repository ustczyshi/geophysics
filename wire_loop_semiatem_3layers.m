%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平电偶源激励，地面或空中观测,层状大地产生的响应分析

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


%%
% [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t);
save('horizontal_electrical_dipole_impulse_shuzhijie','hx_1_impulse','hy_1_impulse','hz_1_impulse');
%% z轴阶跃
figure;
loglog(t.*10^3,u0.*abs(hz_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hz_10),'b:','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% x轴阶跃
figure;
loglog(t.*10^3,u0.*abs(hx_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hx_10),'b:','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
% y 轴阶跃
figure;
loglog(t.*10^3,u0.*abs(hy_01),'r','Linewidth',1);
hold on
loglog(t.*10^3,u0.*abs(hy_10),'b:','Linewidth',1);
grid on;
legend('数值解正阶跃响应','数值解负阶跃响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
%
%% z 轴脉冲
load horizontal_electrical_dipole_impulse_jiexijie.mat;
figure;
plot(t(1:end).*10^3,real(ex_01),'r','Linewidth',2);
hold on
plot(t(1:end).*10^3,(step_ex01),'k:','Linewidth',2);
grid on;
legend('数值解real(ex\_01)','解析解step\_ex01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ex step response'])
xlabel('Time/(ms)')
ylabel('Ex/(V/m)');
figure;
plot(t(1:end).*10^3,real(ey_01),'r','Linewidth',2);
hold on
plot(t(1:end).*10^3,(step_ey01),'k:','Linewidth',2);
grid on;
legend('数值解real(ey\_01)','解析解step\_ey01');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop ey step response'])
xlabel('Time/(ms)')
ylabel('Ey/(V/m)');
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
%{
figure;
plot(t(1:end).*10^3,u0.*(hz_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hz_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hz),'k:','Linewidth',2);
grid on;
legend('数值解hz\_1\_impulse','数值解diff(hz\_01)./dt','解析解impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
plot(t(1:end).*10^3,u0.*(hx_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hx_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hx),'k:','Linewidth',2);
grid on;
legend('数值解hx\_1\_impulse','数值解diff(hx\_01)./dt','解析解impulse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
figure;
plot(t(1:end).*10^3,u0.*(hy_1_impulse),'r','Linewidth',2);
hold on
plot(t(1:end-1).*10^3,u0.*(diff(hy_01)./dt),'b','Linewidth',2);
grid on;
plot(t(1:end).*10^3,u0.*(impulse_hy),'k:','Linewidth',2);
grid on;
legend('数值解hy\_1\_impulse','数值解diff(hy\_01)./dt','解析解impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%}
%%

%% z 轴脉冲
%{
figure;
plot(t(1:end).*10^3,u0.*hz_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,u0.*hz_1_impulse,'b','Linewidth',1);
grid on;
legend('数值解正脉冲响应','数值解负脉冲响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz'])
xlabel('Time/(ms)')
ylabel('Bz/(A/m)');
% y 轴脉冲
figure;
plot(t(1:end).*10^3,u0.*hy_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,-u0.*hy_1_impulse,'b','Linewidth',1);
grid on;
legend('数值解正脉冲响应','数值解负脉冲响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By'])
xlabel('Time/(ms)')
ylabel('By/(A/m)');
% x 轴脉冲
figure;
plot(t(1:end).*10^3,u0.*hx_1_impulse,'r','Linewidth',1);
hold on
plot(t(1:end).*10^3,-u0.*hx_1_impulse,'b','Linewidth',1);
grid on;
legend('数值解正脉冲响应','数值解负脉冲响应');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx'])
xlabel('Time/(ms)')
ylabel('Bx/(A/m)');
%}
%% save data

