function A_Matrix=Get_A_Matrix(N,lambda,gama1,gama2)
%lambda=R2/R1
%gama1=G0/G1
%gama2=G1/G2
A_Matrix=zeros(8*N,8*N);
for i=1:N
    % 1 N hang
    A_Matrix(i,2*i-1)=1;   
    A_Matrix(i,2*i-1+2*N)=-lambda^i;
    A_Matrix(i,2*i-1+4*N)=1;
    % (N+1) 2N hang
    A_Matrix(i+N,2*i)=1;  
    A_Matrix(i+N,2*i+2*N)=-lambda^i;
    A_Matrix(i+N,2*i+4*N)=-1; 
    % (2N+1) 3N hang
    A_Matrix(i+2*N,2*i-1)=1;  
    A_Matrix(i+2*N,2*i-1+2*N)=-gama2*lambda^i;
    A_Matrix(i+2*N,2*i-1+4*N)=-gama2;
    % (3N+1) 4N hang
    A_Matrix(i+3*N,2*i)=1;  
    A_Matrix(i+3*N,2*i+2*N)=-gama2*lambda^i;
    A_Matrix(i+3*N,2*i+4*N)=gama2; 
    % (4N+1) 5N hang
    A_Matrix(i+4*N,2*i-1+2*N)=1;  
    A_Matrix(i+4*N,2*i-1+4*N)=-lambda^i;
    A_Matrix(i+4*N,2*i-1+6*N)=1;
    % (5N+1) 6N hang
    A_Matrix(i+5*N,2*i+2*N)=1;  
    A_Matrix(i+5*N,2*i+4*N)=lambda^i;
    A_Matrix(i+5*N,2*i+6*N)=-1;  
    % (6N+1) 7N hang
    A_Matrix(i+6*N,2*i-1+2*N)=1;  
    A_Matrix(i+6*N,2*i-1+4*N)=lambda^i;
    A_Matrix(i+6*N,2*i-1+6*N)=-gama1;

    A_Matrix(i+7*N,2*i+2*N)=1;  
    A_Matrix(i+7*N,2*i+4*N)=-lambda^i;
    A_Matrix(i+7*N,2*i+6*N)=gama1;     
end

end