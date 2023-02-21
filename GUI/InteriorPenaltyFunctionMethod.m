function InteriorPenaltyFunctionMethod

k=0;
r=1;
c=0.1;
format long;
e=0.001;
N=100;
X0=[2;1];
Xr=[3;1];
disp('0');
disp(X0');

while sqrt(sum((X0-Xr).^2))>e&k<N
    k=k+1;
    syms X1 X2 real;
    f=X1^2+X2^2-r*log(X1-1);
    dfx1=diff(f,X1,1);
    dfx2=diff(f,X2,1);
    dfx1x1=diff(f,X1,2);
    dfx2x2=diff(f,X2,2);
    dfx1x2=diff(dfx1,X2,1);
    grads_x1=inline(vectorize(dfx1),'X1','X2');
    grads_x2=inline(vectorize(dfx2),'X1','X2');
    hession_11=inline(vectorize(dfx1x1),'X1','X2');
    hession_12=inline(vectorize(dfx1x2),'X1','X2');
    hession_22=inline(vectorize(dfx2x2),'X2','X2');
    inv_hession=inv([hession_11(X0(1),X0(2)),hession_12(X0(1),X0(2));hession_12(X0(1),X0(2)),hession_22(X0(1),X0(2))]);
    giads=[grads_x1(X0(1),X0(2));grads_x2(X0(1),X0(2))];
    Xr=X0;
    X0=X0-inv_hession*giads;
    r=r*c;
    disp(k);
    disp(X0');
end
disp('The result is: ');
disp(k);
disp(X0');
    
    


