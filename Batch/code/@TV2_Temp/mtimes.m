function res = mtimes(a,b)

if a.adjoint
    res = adjDz(b);
else
%     res = b(:,:,[2:end,end]) - b;
    
    [rows,cols,timep] = size(b); %'Me'
    res = zeros(rows,cols,timep);
    res(:,:,2:timep-1) = b(:,:,3:timep)-2*b(:,:,2:timep-1)+b(:,:,1:timep-2);
    res(:,:,1) = b(:,:,timep)-2*b(:,:,1)+b(:,:,2);
    res(:,:,timep) = b(:,:,timep-1)-2*b(:,:,timep)+b(:,:,1);
    
end

function y = adjDz(x)
% y= x(:,:,[1,1:end-1]) - x;
% y(:,:,1) = -x(:,:,1);
% y(:,:,end) = x(:,:,end-1);


[rows,cols,timep] = size(x); 
y = zeros(rows,cols,timep);
y(:,:,2:cols-1) = x(:,:,3:timep)-2*x(:,:,2:timep-1)+x(:,:,1:timep-2);
y(:,:,1) = x(:,:,timep)-2*x(:,:,1)+x(:,:,2);
y(:,:,timep) = x(:,:,timep-1)-2*x(:,:,timep)+x(:,:,1);