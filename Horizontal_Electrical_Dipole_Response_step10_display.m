%%  均匀半空间 水平电偶源在地表处的响应分析
clc;clear all;close all;
%% 参数设置
u0 = 4*pi*1e-7;
fs = 1e5;
dt =1./fs;
Ns = 4e3;
Tob = (1:Ns)./fs;
rou = 100;% 电阻率
I = 1;
ds = 1;
m = I.*ds;
xr = 100;
yr = 100;
zr = 0;
position = ['(' num2str(xr) ',' num2str(yr) ',' num2str(zr) ')' ];
%%
[step_ex01,step_ey01,impulse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,xr,yr,zr,Tob);

save('horizontal_electrical_dipole_impulse_jiexijie','step_ex01','step_ey01','impulse_hx','impulse_hy','impulse_hz');
%% 水平电场的阶跃响应
%{
figure;
plot(Tob,step_ex01,'r','linewidth',2);
hold on;
plot(Tob,step_ey01,'b','linewidth',2);
legend('ex','ey');
title(['水平电偶极源在' position '的负阶跃水平电场']);
xlabel('Time/s');
ylabel('E/(V/m)');
%% 磁场的脉冲响应 hx
figure;
plot(Tob,u0.*(impulse_hx),'r','linewidth',2);
% hold on;
% plot(Tob,impulse_hy,'b-.','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title([ '水平电偶极源在' position '产生磁场的脉冲响应Bx']);
xlabel('Time/s');
ylabel('B/(T)');

%%  磁场的脉冲响应 hy
figure;
% plot(Tob,impluse_hx,'r--','linewidth',2);
% hold on;
plot(Tob,u0.*impulse_hy,'b','linewidth',2);
% hold on;
% plot(Tob,impulse_hz,'c','linewidth',2);
% legend('hx','hy','hz');
grid on;
title([ '水平电偶极源在' position '产生磁场的脉冲响应By']);
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
grid on;
title([ '水平电偶极源在' position '产生磁场的脉冲响应Bz']);
xlabel('Time/s');
ylabel('B/(T)');
%}
%% 