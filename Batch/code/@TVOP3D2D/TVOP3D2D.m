function  res = TVOP3D2D()

%res = TVOP3D2D()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res = class(res,'TVOP3D2D');

