function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005
res = cat(4, Dxx(image), Dyy(image), Dxy(image));

function d = Dxx(u)
[rows,cols,sls] = size(u); %'Me'
d = zeros(rows,cols,sls);
d(:,2:cols-1,:) = u(:,3:cols,:)-2*u(:,2:cols-1,:)+u(:,1:cols-2,:);
d(:,1,:) = u(:,cols,:)-2*u(:,1,:)+u(:,2,:);
d(:,cols,:) = u(:,cols-1,:)-2*u(:,cols,:)+u(:,1,:);
return

function d = Dxxt(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(:,2:cols-1,:) = u(:,3:cols,:)-2*u(:,2:cols-1,:)+u(:,1:cols-2,:);
d(:,1,:) = u(:,cols,:)-2*u(:,1,:)+u(:,2,:);
d(:,cols,:) = u(:,cols-1,:)-2*u(:,cols,:)+u(:,1,:);
return

function d = Dyy(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(2:rows-1,:,:) = u(3:rows,:,:)-2*u(2:rows-1,:,:)+u(1:rows-2,:,:);
d(1,:,:) = u(2,:,:)+u(rows,:,:)-2*u(1,:,:);
d(rows,:,:) = u(1,:,:)+u(rows-1,:,:)-2*u(rows,:,:);
return

function d = Dyyt(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(2:rows-1,:,:) = u(3:rows,:,:)-2*u(2:rows-1,:,:)+u(1:rows-2,:,:);
d(1,:,:) = u(2,:,:)+u(rows,:,:)-2*u(1,:,:);
d(rows,:,:) = u(1,:,:)+u(rows-1,:,:)-2*u(rows,:,:);
return

% function d = Dxy(u)
% [rows,cols] = size(u); 
% d = zeros(rows,cols);
% d(1:rows-1,1:cols-1) = u(1:rows-1,1:cols-1)-u(2:rows,1:cols-1)-u(1:rows-1,2:cols)+u(2:rows,2:cols);
% d(rows,1:cols-1) = u(rows,1:cols-1)+u(1,2:cols)-u(1,1:cols-1)-u(rows,2:cols);
% d(1:rows-1,cols) = u(1:rows-1,cols)+u(2:rows,1)-u(1:rows-1,cols)-u(2:rows,1);
% d(rows,cols) = u(rows,cols)+u(1,1)-u(1,cols)-u(rows,1);
% return

function d = Dxy(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(1:rows-1,1:cols-1,:) = u(1:rows-1,1:cols-1,:)-u(2:rows,1:cols-1,:)-u(1:rows-1,2:cols,:)+u(2:rows,2:cols,:);
d(rows,1:cols-1,:) = u(rows,1:cols-1,:)+u(1,2:cols,:)-u(1,1:cols-1,:)-u(rows,2:cols,:);
d(1:rows-1,cols,:) = u(1:rows-1,cols,:)+u(2:rows,1,:)-u(2:rows,cols,:)-u(1:rows-1,1,:);%d(1:rows-1,cols) = u(1:rows-1,cols)+u(2:rows,1)-u(1:rows-1,cols)-u(2:rows,1);
d(rows,cols,:) = u(rows,cols,:)+u(1,1,:)-u(1,cols,:)-u(rows,1,:);
return

function d = Dxyt(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(2:rows,2:cols,:) = u(2:rows,2:cols,:)-u(2:rows,1:cols-1,:)-u(1:rows-1,2:cols,:)+u(1:rows-1,1:cols-1,:);
d(1,2:cols,:) = u(1,2:cols,:)+u(rows,1:cols-1,:)-u(1,1:cols-1,:)-u(rows,2:cols,:);
d(2:rows,1,:) = u(2:rows,1,:)+u(1:rows-1,cols,:)-u(1:rows-1,1,:)-u(2:rows,cols,:);
d(1,1,:) = u(1,1,:)+u(rows,cols,:)-u(1,cols,:)-u(rows,1,:);
return

function d = Dyx(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(1:rows-1,1:cols-1,:) = u(1:rows-1,1:cols-1,:)-u(2:rows,1:cols-1,:)-u(1:rows-1,2:cols,:)+u(2:rows,2:cols,:);
d(rows,1:cols-1,:) = u(rows,1:cols-1,:)+u(1,2:cols,:)-u(1,1:cols-1,:)-u(rows,2:cols,:);
d(1:rows-1,cols,:) = u(1:rows-1,cols,:)+u(2:rows,1,:)-u(1:rows-1,cols,:)-u(2:rows,1,:);
d(rows,cols,:) = u(rows,cols,:)+u(1,1,:)-u(1,cols,:)-u(rows,1,:);
return

function d = Dyxt(u)
[rows,cols, sls] = size(u); 
d = zeros(rows,cols, sls);
d(2:rows,2:cols,:) = u(2:rows,2:cols,:)-u(2:rows,1:cols-1,:)-u(1:rows-1,2:cols,:)+u(1:rows-1,1:cols-1,:);
d(1,2:cols,:) = u(1,2:cols,:)+u(rows,1:cols-1,:)-u(1,1:cols-1,:)-u(rows,2:cols,:);
d(2:rows,1,:) = u(2:rows,1,:)+u(1:rows-1,cols,:)-u(1:rows-1,1,:)-u(2:rows,cols,:);
d(1,1,:) = u(1,1,:)+u(rows,cols,:)-u(1,cols,:)-u(rows,1,:);
return



