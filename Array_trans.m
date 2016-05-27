%设置取样点lambda（i）。并将lambda和频率扩展为二维矩阵，
%使不同行频率不同，同一行频率相同；不同列lambda不同，同一列lambda相同。
% lambda: 行向量；
% freq：行向量；
function [lmd_array,fre_array] = Array_trans(lambda, freq)
if (size(lambda,1)>1) &&(size(lambda,2)>1)
    error('input lambda is not a vector');  % 输入非向量则报错
else
    if (size(lambda,2)==1) %% 如果输入为列向量，转化为行向量
        lambda = lambda';
    end
end
    
if (size(freq,1)>1) &&(size(freq,2)>1)
    error('input freq is not a vector'); 
else
    if (size(freq,2)==1)
        freq = freq';
    end
end
col= length(lambda);
row = length(freq);
lmd_array = repmat(lambda,row,1);%不同列lambda不同，同一列lambda相同。
fre_array = repmat(freq',1,col);%使不同行频率不同，同一行频率相同；
% lmd_array = ones(row,1) * lambda;
% fre_array = freq' * ones(1,col);
end