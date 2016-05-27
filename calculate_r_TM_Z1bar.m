function rtm_z1bar=calculate_r_TM_Z1bar(lambda,freq)
%% 计算TM波在分层大地中的反射系数,用于电性源分反演,对于大多数勘查地球物理方法rtm=1
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
kk = -1i * 2 * pi * frequency_Array * mu_0;%使不同行频率不同，同一行频率相同
% u1_star =  (lambda_2 - kk .* sigma(N)).^0.5;% 获得最底层的u_n
%%   新添加项-------------------------------------20160511
z1_star =  (lambda_2 - kk .* sigma(N)).^0.5./sigma(N);
for k= N-1:-1:1 %% n层只需要递推n-1次
     u_n = (lambda_2 - kk.*sigma(k)).^0.5;% 第n层的u_n；
    z_n = u_n./sigma(k);% 第n层的z_n,新添加项------------20160511
%     u1_star = u_n .*( u1_star + u_n .* exp(-2*u_n * d(k)))./(u_n + u1_star .* exp(-2*u_n * d(k)));
    z1_star = z_n .*( z1_star + z_n .* exp(-2*u_n * d(k)))./(z_n + z1_star .* exp(-2*u_n * d(k)));
end
rtm_z1bar = z1_star;

end
