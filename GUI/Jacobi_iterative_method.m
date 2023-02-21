function Jacobi_iterative_method

format long;
A=input('please input coefficient matrix row by row: ');
b=input('please input the right value of the function by column: ');
X0=input('please input the initial value of vector X by column: ');
e=input('please input permitted error limited value: ');
N=input('please input the MAX iterative times: ');

K=0;
n=length(A(1,:));
D=diag(diag(A));
B=eye(n)-D^-1*A;
g=D^-1*b;
X1=B*X0+g;
K=K+1;
fprintf('\n');
fprintf('The initial value of X as follows: ');
fprintf(1,'%f  ', X0');
fprintf('\n');
fprintf(1,'The %d times iterative X value is:',K);
fprintf(1,'%15.12f',X1);
fprintf('\n');

while sqrt(sum((X1-X0).^2))>e&K<N
    K=K+1;
    X0=X1;
    X1=B*X0+g;
    fprintf(1,'The %d times iterative X value is: ',K);
    fprintf(1,'%15.12f',X1);
    fprintf('\n');
end

if K<N
    fprintf('\n');
    disp(['The result value of X is: ' num2str(X1')]);
    disp(['The iterative times is: ' num2str(K)]);
    fprintf('\n');
else
    fprintf('\n');
    disp('Warning: It is not a convergence function!');
    fprintf('\n');
end

