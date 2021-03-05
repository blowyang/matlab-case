function A_Matrix=Get_A_Matrix(N,K,R1,R2,G0,G2)
%N 级数截断次数
%K 界面相分层数
h=(R1-R2)/K;
A_Matrix=zeros((2+2*K)*2*N,(2+2*K)*2*N);
for j=1:N
    lambda=R2/(R2+h);
    Rk=R2+0.5*h;
    Gk=Get_G1(G0,G2,R1,R2,Rk);
    gama2=Gk/G2;
    % 1 N hang
    A_Matrix(j,2*j-1)=1;   
    A_Matrix(j,2*j-1+2*N)=-lambda^j;
    A_Matrix(j,2*j-1+4*N)=1;
    % (N+1) 2N hang
    A_Matrix(j+N,2*j)=1;  
    A_Matrix(j+N,2*j+2*N)=-lambda^j;
    A_Matrix(j+N,2*j+4*N)=-1; 
    % (2N+1) 3N hang
    A_Matrix(j+2*N,2*j-1)=1;  
    A_Matrix(j+2*N,2*j-1+2*N)=-gama2*lambda^j;
    A_Matrix(j+2*N,2*j-1+4*N)=-gama2;
    % (3N+1) 4N hang
    A_Matrix(j+3*N,2*j)=1;  
    A_Matrix(j+3*N,2*j+2*N)=-gama2*lambda^j;
    A_Matrix(j+3*N,2*j+4*N)=gama2; 
    for k=1:K-1
        A_Matrix(j+4*k*N,2*j-1+2*N+4*N*(k-1))=1;
        A_Matrix(j+4*k*N,2*j-1+2*N+4*N*(k-1)+2*N)=-((R2+(k-1)*h)/(R2+k*h))^j;
        A_Matrix(j+4*k*N,2*j-1+2*N+4*N*(k-1)+2*N+2*N)=-((R2+k*h)/(R2+(k+1)*h))^j;
        A_Matrix(j+4*k*N,2*j-1+2*N+4*N*(k-1)+2*N+2*N+2*N)=1;
        
        A_Matrix(j+4*k*N+N,2*j+2*N+4*N*(k-1))=1;
        A_Matrix(j+4*k*N+N,2*j+2*N+4*N*(k-1)+2*N)=((R2+(k-1)*h)/(R2+k*h))^j;
        A_Matrix(j+4*k*N+N,2*j+2*N+4*N*(k-1)+2*N+2*N)=-((R2+k*h)/(R2+(k+1)*h))^j;
        A_Matrix(j+4*k*N+N,2*j+2*N+4*N*(k-1)+2*N+2*N+2*N)=-1; 
        
        Rk=R2+(k-1)*h+0.5*h;
        Rk_=R2+k*h+0.5*h;        
        gamak=Get_G1(G0,G2,R1,R2,Rk_)/Get_G1(G0,G2,R1,R2,Rk);
        A_Matrix(j+4*k*N+2*N,2*j-1+2*N+4*N*(k-1))=1;
        A_Matrix(j+4*k*N+2*N,2*j-1+2*N+4*N*(k-1)+2*N)=((R2+(k-1)*h)/(R2+k*h))^j;
        A_Matrix(j+4*k*N+2*N,2*j-1+2*N+4*N*(k-1)+2*N+2*N)=-gamak*((R2+k*h)/(R2+(k+1)*h))^j;
        A_Matrix(j+4*k*N+2*N,2*j-1+2*N+4*N*(k-1)+2*N+2*N+2*N)=-gamak;  
        
        A_Matrix(j+4*k*N+3*N,2*j+2*N+4*N*(k-1))=1;
        A_Matrix(j+4*k*N+3*N,2*j+2*N+4*N*(k-1)+2*N)=-((R2+(k-1)*h)/(R2+k*h))^j;
        A_Matrix(j+4*k*N+3*N,2*j+2*N+4*N*(k-1)+2*N+2*N)=-gamak*((R2+k*h)/(R2+(k+1)*h))^j;
        A_Matrix(j+4*k*N+3*N,2*j+2*N+4*N*(k-1)+2*N+2*N+2*N)=gamak;
    end
    % (4N+1) 5N hang
    lambda=(R1-h)/R1;
    Rk=R1-0.5*h;
    gama1=G0/Get_G1(G0,G2,R1,R2,Rk);
    A_Matrix(j+4*K*N,2*j-1+4*K*N-2*N)=1;  
    A_Matrix(j+4*K*N,2*j-1+4*K*N)=-lambda^j;
    A_Matrix(j+4*K*N,2*j-1+4*K*N+2*N)=1;
    % (5N+1) 6N hang
    A_Matrix(j+4*K*N+N,2*j+4*K*N-2*N)=1;
    A_Matrix(j+4*K*N+N,2*j+4*K*N)=lambda^j;
    A_Matrix(j+4*K*N+N,2*j+4*K*N+2*N)=-1;  
    % (6N+1) 7N hang
    A_Matrix(j+4*K*N+2*N,2*j-1+4*K*N-2*N)=1;  
    A_Matrix(j+4*K*N+2*N,2*j-1+4*K*N)=lambda^j;
    A_Matrix(j+4*K*N+2*N,2*j-1+4*K*N+2*N)=-gama1;

    A_Matrix(j+4*K*N+3*N,2*j+4*K*N-2*N)=1;  
    A_Matrix(j+4*K*N+3*N,2*j+4*K*N)=-lambda^j;
    A_Matrix(j+4*K*N+3*N,2*j+4*K*N+2*N)=gama1;     
end

end