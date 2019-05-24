function [ Accuarcy] = EvaluationFcn( C_hat,outer_mat)
%EVALUATIONFCN 指标计算
%   首先将向量化为对角线为零，再求解TPR和准确率；
SIZE = size(outer_mat,1);
[m,n]= size(outer_mat);
epsi = 0.1;                                     %阈值
for i = 1:SIZE
    C_hat(i,i) = 1;
end

adjRecon = reshape(C_hat,SIZE * SIZE,1);          %重构矩阵换成向量形式
adjReal  = reshape(outer_mat,SIZE*SIZE,1);      %真实矩阵换成向量形式

for i = 1:size(adjRecon)
    if ( abs(adjRecon(i)-1 ) <= epsi)
        adjRecon(i) = 1;
    elseif ( abs(adjRecon(i)-0 ) <= epsi)
            adjRecon(i) = 0;
    end
end

% 计算Accuarcy
acc = 0;
tp = 0;
for i = 1:size(adjReal)
    if (adjReal(i) == adjRecon(i))
        acc = acc +1;
    end
end
tp  = size(adjReal,1);
Accuarcy = acc/tp;

% for i = 1:size(C_hat,1)
%     for j = 1:size(C_hat,1)
%         x(i,j) = j/100 -0.5;
%     end
% end
% x1 = reshape(x,SIZE*SIZE,1);
% % y1 = reshape(adjRecon,SIZE*SIZE,1);
% scatter(x1,adjRecon);

end