function y=Fast_Hankel(r,sum,J_zero)
%r为收发距，
load parameters.txt;% 读取参数
m = parameters(1,1); % 读取磁矩
y = m/r  *  sum * J_zero';% 用矩阵乘法计算数值积分，y是个列向量，每一列对应不同的频率
%./r是由汉克尔变换过程中的转换构造引出的
end