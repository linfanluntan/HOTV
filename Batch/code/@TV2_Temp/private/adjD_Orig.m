function res = adjD(y)

res = zeros(size(y,1),size(y,2));

%y1 = ones(imsize)*y(1)/sqrt(prod(imsize));
%yx = (reshape(y(2:prod(imsize)+1), imsize(1), imsize(2)));
%yy = (reshape(y(prod(imsize)+2:end), imsize(1), imsize(2)));

res = adjDx(y(:,:,1)) + adjDy(y(:,:,2)) + adjDd(y(:,:,3));

end


function res = adjDy(x)

res1 = x(:,[1,1:end-1]) - x;
res1(:,1) = -x(:,1);
res1(:,end) = x(:,end-1);

%
res2 = x(:,[2:end end]) - x;
res2(:,1) = x(:,2);
res2(:,end) = -x(:,end);

%
res = res1 + res2;

end

function res = adjDx(x)

res1 = x([1,1:end-1],:) - x;
res1(1,:) = -x(1,:);
res1(end,:) = x(end-1,:);

%
res2 = x([2:end end],:) - x;
res2(1,:) = x(2,:);
res2(end,:) = -x(end,:);

%
res = res1 + res2;

end


function res = adjDd(x)

res1 = x([2:end end],:) - x;
res1(1,:) = x(2,:);
res1(end,:) = -x(end,:);

res2 = x([2:end end],[2:end end]) - x(:,[2:end end]);
res2(1,[2:end end]) = x(2,[2:end end]);
res2(end,[2:end end]) = -x(end,[2:end end]);
res2(:,end) = 0;

res = res1 - res2;

end

