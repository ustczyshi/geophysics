function rte=calculate_r_TE(lambda,freq)
%% 计算TE波在分层大地中的反射系数
%定义参数
mu_0 = 4*pi*1e-7;
%读取参数
load parameters.txt;
sigma = parameters(:,2)';% 电导率 1~n各元素代表1~n层的电导率 列向量转化为行向量
d =  parameters(:,3)'; % 层厚
N=length(sigma);
[lambda_Array,frequency_Array] = Array_trans(lambda,freq); % 将向量转化为矩阵
% for循环中多次计算部分单独提出，保存中间数据，减小运算量
lambda_2 = lambda_Array.^2; 
% kk = -1i * 2 * pi * frequency_Array * mu_0;
% u1_star =  (lambda_2 - kk * 1./rho(N)).^0.5;% 
kk = -1i * 2 * pi * frequency_Array * mu_0;
u1_star =  (lambda_2 - kk .* sigma(N)).^0.5;% 
for k= N-1:-1:1 %% n层只需要递推n-1次
    u_n = (lambda_2 - kk.*sigma(k)).^0.5;% 第n层的u_n；
    u1_star = u_n .*( u1_star + u_n .* exp(-2*u_n * d(k)))./(u_n + u1_star .* exp(-2*u_n * d(k)));
end
rte = (lambda_Array-u1_star)./(lambda_Array+u1_star);

end
