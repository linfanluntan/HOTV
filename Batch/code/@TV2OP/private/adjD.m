function res = adjD(y)
%'Div2, I run'
%res = zeros(size(y,1),size(y,2));

%y1 = ones(imsize)*y(1)/sqrt(prod(imsize));
%yx = (reshape(y(2:prod(imsize)+1), imsize(1), imsize(2)));
%yy = (reshape(y(prod(imsize)+2:end), imsize(1), imsize(2)));

%res = adjDx(y(:,:,1)) + adjDy(y(:,:,2)) + adjDd(y(:,:,3));
res = +Dxxt(y(:,:,1)) + Dyyt(y(:,:,2)) + Dxyt(y(:,:,3)) ;


function d = Dxx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols-1) = u(:,3:cols)-2*u(:,2:cols-1)+u(:,1:cols-2);
d(:,1) = u(:,cols)-2*u(:,1)+u(:,2);
d(:,cols) = u(:,cols-1)-2*u(:,cols)+u(:,1);
return

function d = Dxxt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols-1) = u(:,3:cols)-2*u(:,2:cols-1)+u(:,1:cols-2);
d(:,1) = u(:,cols)-2*u(:,1)+u(:,2);
d(:,cols) = u(:,cols-1)-2*u(:,cols)+u(:,1);
return

function d = Dyy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows-1,:) = u(3:rows,:)-2*u(2:rows-1,:)+u(1:rows-2,:);
d(1,:) = u(2,:)+u(rows,:)-2*u(1,:);
d(rows,:) = u(1,:)+u(rows-1,:)-2*u(rows,:);
return

function d = Dyyt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows-1,:) = u(3:rows,:)-2*u(2:rows-1,:)+u(1:rows-2,:);
d(1,:) = u(2,:)+u(rows,:)-2*u(1,:);
d(rows,:) = u(1,:)+u(rows-1,:)-2*u(rows,:);
return

function d = Dxy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(1:rows-1,1:cols-1) = u(1:rows-1,1:cols-1)-u(2:rows,1:cols-1)-u(1:rows-1,2:cols)+u(2:rows,2:cols);
d(rows,1:cols-1) = u(rows,1:cols-1)+u(1,2:cols)-u(1,1:cols-1)-u(rows,2:cols);
d(1:rows-1,cols) = u(1:rows-1,cols)+u(2:rows,1)-u(1:rows-1,cols)-u(2:rows,1);
d(rows,cols) = u(rows,cols)+u(1,1)-u(1,cols)-u(rows,1);
return

function d = Dxyt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,2:cols) = u(2:rows,2:cols)-u(2:rows,1:cols-1)-u(1:rows-1,2:cols)+u(1:rows-1,1:cols-1);
d(1,2:cols) = u(1,2:cols)+u(rows,1:cols-1)-u(1,1:cols-1)-u(rows,2:cols);
d(2:rows,1) = u(2:rows,1)+u(1:rows-1,cols)-u(1:rows-1,1)-u(2:rows,cols);
d(1,1) = u(1,1)+u(rows,cols)-u(1,cols)-u(rows,1);
return

function d = Dyx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(1:rows-1,1:cols-1) = u(1:rows-1,1:cols-1)-u(2:rows,1:cols-1)-u(1:rows-1,2:cols)+u(2:rows,2:cols);
d(rows,1:cols-1) = u(rows,1:cols-1)+u(1,2:cols)-u(1,1:cols-1)-u(rows,2:cols);
d(1:rows-1,cols) = u(1:rows-1,cols)+u(2:rows,1)-u(1:rows-1,cols)-u(2:rows,1);
d(rows,cols) = u(rows,cols)+u(1,1)-u(1,cols)-u(rows,1);
return

function d = Dyxt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,2:cols) = u(2:rows,2:cols)-u(2:rows,1:cols-1)-u(1:rows-1,2:cols)+u(1:rows-1,1:cols-1);
d(1,2:cols) = u(1,2:cols)+u(rows,1:cols-1)-u(1,1:cols-1)-u(rows,2:cols);
d(2:rows,1) = u(2:rows,1)+u(1:rows-1,cols)-u(1:rows-1,1)-u(2:rows,cols);
d(1,1) = u(1,1)+u(rows,cols)-u(1,cols)-u(rows,1);
return

