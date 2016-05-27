%%  均匀半空间 水平电偶源在地表处的响应分析
clc;clear all;close all;
%% 参数设置
u0 = 4*pi*1e-7;
fs = 1e7;
dt =1./fs;
Ns = 2e4;
Tob = (1:Ns)./fs;
rou = 100;% 电阻率
I = 1;
ds = 1;
m = I.*ds;
xr = 100;
yr = 0;
zr = 0;
%%
[step_ex10,step_ey10,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,xr,yr,zr,Tob);
%% 水平电场的阶跃响应
figure;
plot(Tob,step_ex10,'r','linewidth',2);
hold on;
plot(Tob,step_ey10,'b','linewidth',2);
legend('ex','ey');
title('水平电偶极源在(0,100,0)的阶跃水平电场');
xlabel('Time/s');
ylabel('E/(V/m)');
%% 磁场的脉冲响应 hx
figure;
plot(Tob,u0.*impulse_hx,'r','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
title('水平电偶极源在(0,100,0)产生磁场的脉冲响应Bx');
xlabel('Time/s');
ylabel('B/(T)');
% 磁场的阶跃响应 hx
%%  磁场的脉冲响应 hy
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
plot(Tob,u0.*impulse_hy,'b','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title('水平电偶极源在(0,100,0)产生磁场的脉冲响应By');
xlabel('Time/s');
ylabel('B/(T)');
%% 磁场的脉冲响应 hz
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
plot(Tob,impulse_hz.*u0,'r','linewidth',2);
% legend('hx','hy','hz');
title('水平电偶极源在(0,100,0)产生磁场的脉冲响应Bz');
xlabel('Time/s');
ylabel('B/(T)');

%% 