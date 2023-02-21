function [x] = TPTVTV2TV2(imdata_org,param)

par1=1;   % the Lipschitz of gradient of f as indicated in the doc as Beta 2017-4-11
par2=.5;    % the Rou constant as in the doc, 0<Rou<1, 2017-4-11

par3=4;%4+8+48;%8;    % the estimation of operators norm indicated in line H as in the doc,||sum(L*L)||   2017-4-11

par4 = 1;   % the sigma constant as in the doc 2017-4-11
par40 = 1;   % the sigma constant as in the doc 2017-4-11
par4a = 1;     % the Lipschitz of gradient of f indicated in line H as in the doc 2017-4-11
par4b = 1;

% par5=.25/((par3)*par4+par1/2);% .5/((par3)*par4+par1/2);%.9/((par3)*par4+par1/2);
par5=.001/((par3)*par4+par1/2);    % par 5 is the Tao constant as in the doc, however, since we can not estimate the operator norm well, we need to adjust numerator in the right to make iteration convergent   2017-4-11
                              % the recommanded numerator value by testing for cardiac CINE is 0.25 as we have reported   2017-4-11
                               % the bigger the numerator value the fast the convergence as we have reported   2017-4-11
                               % however if the numerator value is bigger then 0.3 the iteration may failed as we have reported   2017-4-11
                               
par5=.01/((par3)*par4+par1/2); 
par5=.02/((par3)*par4+par1/2); 
par5=.04/((par3)*par4+par1/2); 
par5=.08/((par3)*par4+par1/2); 
par5=.16/((par3)*par4+par1/2);
par5=.25/((par3)*par4+par1/2);
par5=.3/((par3)*par4+par1/2);
%  % below is for the definition of projector which is the proximate operator of indicate function   2017-4-11
l=-Inf;% the feasible range of variables to be optimarized is -inf to +inf  2017-4-11
u=Inf;
if((l==-Inf)&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&(l==-Inf))
    project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&isfinite(l))&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end

OptA=project; % prox operator in line C as indicated in the doc  2017-4-11
OptB = @(x)param.TV_dim1*x;  % temporal TV operator as indicated as one of Lm in the doc  2017-4-11
OptC = @(x)param.TV_dim1'*x; % adjoint of temporal TV operator, which need to be used as one of Lm* in the doc  2017-4-11
OptB0 = @(x)param.TV_dim1*x;  % temporal TV operator as indicated as one of Lm in the doc  2017-4-11
OptC0 = @(x)param.TV_dim1'*x; % adjoint of temporal TV operator, which need to be used as one of Lm* in the doc  2017-4-11

OptB1 = @(x)param.TVOP3D2D*x;% spatial TV operator as indicated as one of Lm in the doc  2017-4-11
OptC1 = @(x)param.TVOP3D2D'*x;% adjoint of spatial TV operator, which need to be used as one of Lm* in the doc  2017-4-11

OptB2 = @(x)param.TV2OP3D2D*x;% spatial TV2 operator as indicated as one of Lm in the doc  2017-4-11
OptC2 = @(x)param.TV2OP3D2D'*x;% adjoint of spatial TV2 operator, which need to be used as one of Lm* in the doc  2017-4-11


    y0 = imdata_org;

    par6 = param.TVWeight_dim1;% lamda 1 for temporal TV  2017-4-11
    par60 = param.TV2Weight_dim1;% lamda 1 for temporal TV  2017-4-11
    par6a = param.TVOPWeight;% lamda 2 for spatial TV  2017-4-11
    par6b = param.TV2OPWeight;% lamda 3 for spatial TV2  2017-4-11
    
    OptD =  @(u,par4) par6*u./max(par6, abs(u));% prox operator as indicated in line B for temporal TV  2017-4-11
    OptD0 =  @(u,par4) par6*u./max(par60, abs(u));% prox operator as indicated in line B for temporal TV  2017-4-11
    OptD1 =  @(u,par4) par6a*u./max(par6a, abs(u));% prox operator as indicated in line B for spatial TV  2017-4-11
    OptD2 =  @(u,par4) par6b*u./max(par6b, abs(u));% prox operator as indicated in line B for spatial TV2  2017-4-11
    t=1;
    
    x=y0;  % all initialization   2017-4-11
    u1=OptB(y0)*0;
    u10=OptB0(y0)*0;
    u2=OptB1(y0)*0;
    u3=OptB2(y0)*0;
    con_curve=[];

    for I_0 = 1:param.nite
%         if I_0 >=6
%             par5=.2/((par3)*par4+par1/2);
%         end
%         if I_0 >=13
%             par5=.3/((par3)*par4+par1/2);
%         end
%         if I_0 >=13+(26-15)
%             par5=.2/((par3)*par4+par1/2);
%         end
            
        
        u1p=OptD(u1+par4*OptB(x),par4); % complete pt 1 of line B as in the doc  2017-4-11
        u1old=u1;
        u1=par2*u1p+(1-par2)*u1;% complete pt 2 of line B as in the doc  2017-4-11
        
        u1p0=OptD0(u10+par4*OptB0(x),par40); % complete pt 1 of line B as in the doc  2017-4-11
        u1old0=u10;
        u10=par2*u1p0+(1-par2)*u10;% complete pt 2 of line B as in the doc  2017-4-11
        
        u2p=OptD1(u2+par4a*OptB1(x),par4a);% complete pt 1 of line B as in the doc  2017-4-11
        u2old=u2;
        u2=par2*u2p+(1-par2)*u2;% complete pt 2 of line B as in the doc  2017-4-11
        
        u3p=OptD2(u3+par4b*OptB2(x),par4b);% complete pt 1 of line B as in the doc  2017-4-11
        u3old=u3;
        u3=par2*u3p+(1-par2)*u3;% complete pt 2 of line B as in the doc  2017-4-11
        
        
%          xp=ProxG(x-tau*(ifft2c(conj(R).*(R.*fft2c(x)-yk0)))-tau*(KS1(2*u1p-u1old)+KS2(2*u2p-u2old)));%tau is not included as an indenpendt variable
         xp=OptA(x-par5*((x-param.y))-par5*(OptC(2*u1p-u1old)+OptC0(2*u1p0-u1old0)+OptC1(2*u2p-u2old)+OptC2(2*u3p-u3old)));% complete line C as in the doc  2017-4-11
        xold=x;
        x=par2*xp+(1-par2)*x;% complete line D as in the doc  2017-4-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% below should be uncommented if you need FISTA acc , they are indicated in line E,F,G as in the doc    2017-4-11    
        par7 = (1+sqrt(1+4*t^2))/2;
        par8=(t-1)/(par7);
        t = par7;
        
        x=x+par8*(x-xold);
        u1=u1+par8*(u1-u1old);
        u10=u10+par8*(u10-u1old0);
        u3=u3+par8*(u3-u3old);
        u2=u2+par8*(u2-u2old);
        
%         if mod(I_0,5)==0
%             fname=['swap_int_iterations' '_' num2str(I_0) '.mat'];
%             save(fname,'x');
%         end
        con_curve=[con_curve objective(x,param)];  % below should be commented, it's purly for save temp results      2017-4-11 
        if I_0 == param.nite
            save('it_test001','con_curve');
        end
        

    end
    
return


function res = objective(x,param) %********************************** this can be totally ignored, purely for getting the cost function changes with the iterations        2017-4-11 

% L2-norm part
w=x-param.y;
L2Obj=w(:)'*w(:);

% TV part along time
if param.TVWeight_dim1
    w = param.TV_dim1*(x); 
    TVObj_dim1 = sum((w(:).*conj(w(:))).^(1/2));
else
    TVObj_dim1 = 0;
end

if param.TV2Weight_dim1
    w = param.TV2_dim1*(x); 
    TVObj_dim10 = sum((w(:).*conj(w(:))).^(1/2));
else
    TVObj_dim10 = 0;
end

% TV part along respiration
if param.TVOPWeight
    w = param.TVOP3D2D*(x); 
    TVObj_dim2 = sum((w(:).*conj(w(:))).^(1/2));
else
    TVObj_dim2 = 0;
end
if param.TV2OPWeight
    w = param.TV2OP3D2D*(x); 
    TVObj_dim3 = sum((w(:).*conj(w(:))).^(1/2));
else
    TVObj_dim3 = 0;
end

res1=L2Obj+param.TVWeight_dim1*TVObj_dim1+param.TV2Weight_dim1*TVObj_dim10+param.TVOPWeight*TVObj_dim2+param.TV2OPWeight*TVObj_dim3;
res0=0.5*L2Obj+param.TVWeight_dim1*TVObj_dim1+param.TV2Weight_dim1*TVObj_dim10+param.TVOPWeight*TVObj_dim2+param.TV2OPWeight*TVObj_dim3;
res=[res0 res1];


% function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)
% 
% % precalculates transforms to make line search cheap
% 
% FTXFMtx = params.FT*(params.XFM'*x);
% FTXFMtdx = params.FT*(params.XFM'*dx);
% 
% if params.TVWeight
%     DXFMtx = params.TV*(params.XFM'*x);
%     DXFMtdx = params.TV*(params.XFM'*dx);
% else
%     DXFMtx = 0;
%     DXFMtdx = 0;
% end
% 
% 
% 
% 
% 
% function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params);
% %calculated the objective function
% 
% p = params.pNorm;
% 
% obj = FTXFMtx + t*FTXFMtdx - params.data;
% obj = obj(:)'*obj(:);
% 
% if params.TVWeight
%     w = DXFMtx(:) + t*DXFMtdx(:);
%     TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
% else
%     TV = 0;
% end
% 
% if params.xfmWeight
%    w = x(:) + t*dx(:); 
%    XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
% else
%     XFM=0;
% end
% 
% 
% 
% TV = sum(TV.*params.TVWeight(:));
% XFM = sum(XFM.*params.xfmWeight(:));
% RMS = sqrt(obj/sum(abs(params.data(:))>0));
% 
% res = obj + (TV) + (XFM) ;

