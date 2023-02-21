function LinearCofficientMethod
X0=[-1;-1;-1];
format long;
M0=1;
C=8;
e1=0.001;
e2=0.001;
k=0;
t=0;
kk=0;
k1=0;
N=100;
N1=50;
n=length(X0);
while t<N1
    syms X1 X2 X3 real;
    gux1=X1;
    gux2=X2;
    gux3=X3;
    gux4=X1+X2+X3-6;
    gux1_func=inline(vectorize(gux1),'X1','X2','X3');
    gux2_func=inline(vectorize(gux2),'X1','X2','X3');
    gux3_func=inline(vectorize(gux3),'X1','X2','X3');
    if gux1_func(X0(1),X0(2),X0(3))>=0
        gux1=0;
    end
    if gux2_func(X0(1),X0(2),X0(3))>=0
        gux2=0;
    end
    if gux3_func(X0(1),X0(2),X0(3))>=0
        gux3=0;
    end
    fx=0.36*((X1-1)^2+(X2-2)^2+(X3-3)^2)+0.64*(X1^2+2*X2^2+3*X3^2)+M0*(gux1^2+gux2^2+gux3^2+gux4^2);
    fx_func=inline(vectorize(fx),'X1','X2','X3');
    Xr=X0;
    t=t+1;
    while kk<N
        kk=kk+1;
        if k==0
            H0=eye(3);
            g0_1=diff(fx,'X1',1);
            g0_2=diff(fx,'X2',1);
            g0_3=diff(fx,'X3',1);
            g0_1_func=inline(vectorize(g0_1),'X1','X2','X3');
            g0_2_func=inline(vectorize(g0_2),'X1','X2','X3');
            g0_3_func=inline(vectorize(g0_3),'X1','X2','X3');
            g0=[g0_1_func(X0(1),X0(2),X0(3));g0_2_func(X0(1),X0(2),X0(3));g0_3_func(X0(1),X0(2),X0(3))];
        end
        Xi=X0;
        Sk=H0*g0;
        syms r real;
        X0=X0+r*Sk;
        fx_r=vpa(fx_func(X0(1),X0(2),X0(3)),5);
        fx_r_func=diff(fx_r,r,1);
        fx_r_func_func=inline(vectorize(fx_r_func),'r');
        options=optimset('Display','off');
        r=fzero(fx_r_func_func,1,options);
        X0=Xi+r*Sk;
        g0_1=diff(fx,'X1',1);
        g0_2=diff(fx,'X2',1);
        g0_3=diff(fx,'X3',1);
        g0_1_func=inline(vectorize(g0_1),'X1','X2','X3');
        g0_2_func=inline(vectorize(g0_2),'X1','X2','X3');
        g0_3_func=inline(vectorize(g0_3),'X1','X2','X3');
        g0=[g0_1_func(X0(1),X0(2),X0(3));g0_2_func(X0(1),X0(2),X0(3));g0_3_func(X0(1),X0(2),X0(3))];
        if sqrt(sum(g0.^2))<=e1
            break;
        end
        if k==n
            k=0;
        else
            g1=[g0_1_func(Xi(1),Xi(2),Xi(3));g0_2_func(Xi(1),Xi(2),Xi(3));g0_3_func(Xi(1),Xi(2),Xi(3))];
            delta_g=g0-g1;
            delta_x=X0-Xi;
            Ak=delta_x*(delta_x')/(delta_x'*delta_g);
            Bk=H0*delta_g*(H0*delta_g)'*inv(delta_g'*H0*delta_g);
            H0=H0+Ak-Bk;
            k=k+1;
        end
    end
    if sqrt(sum((Xi-X0).^2))<e1&abs((fx_func(Xi(1),Xi(2),Xi(3))-fx_func(X0(1),X0(2),X0(3)))/fx_func(Xi(1),Xi(2),Xi(3)))<e2
        break;
    else
        M0=M0*C;
    end
end
disp(t);
disp(X0);  