% function y = GS_Trans2(t,h_v)
% G_S=load ('G_S.txt')';% 读取GS变换系数  G_S为行向量
function y = GS_Trans2(t,h_v,G_S)
y = log(2)./t .* G_S * h_v;%  脉冲响应
end