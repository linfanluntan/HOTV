function  res = TV2OP3D3D()

%res = TV2OP3D2D()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res = class(res,'TV2OP3D3D');

