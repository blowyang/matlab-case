function A_Matrix=Get_A_Matrix(N,R1,R2,x1,x2,G0,G1,G2)
A_Matrix=zeros(8*N,8*N);
gama10=G1/G0;
gama20=G2/G0;
for j=1:N
    % 1 N hang
    A_Matrix(j,2*j-1)=1;   
    A_Matrix(j,2*j-1+4*N)=1;
%     for k=1:N
%         A_Matrix(j,2*k-1+6*N)=-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j;
%     end
    for k=1:N
        A_Matrix(j,2*k-1+6*N)=real(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
        A_Matrix(j,2*k+6*N)=-imag(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
    end
    % (N+1) 2N hang
    A_Matrix(j+N,2*j)=1;  
    A_Matrix(j+N,2*j+4*N)=-1;
%     for k=1:N
%         A_Matrix(j+N,2*k+6*N)=-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j;
%     end
    for k=1:N
        A_Matrix(j+N,2*k-1+6*N)=imag(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
        A_Matrix(j+N,2*k+6*N)=real(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
    end
    % (2N+1) 3N hang
    
    A_Matrix(j+2*N,2*j-1)=gama10;  
    A_Matrix(j+2*N,2*j-1+4*N)=-1;
%     for k=1:N
%         A_Matrix(j+2*N,2*k-1+6*N)=-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j;
%     end
    for k=1:N
        A_Matrix(j+2*N,2*k-1+6*N)=real(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
        A_Matrix(j+2*N,2*k+6*N)=-imag(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
    end
    % (3N+1) 4N hang
    A_Matrix(j+3*N,2*j)=gama10;  
    A_Matrix(j+3*N,2*j+4*N)=1;
%     for k=1:N
%         A_Matrix(j+3*N,2*k+6*N)=-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j;
%     end
    for k=1:N
        A_Matrix(j+3*N,2*k-1+6*N)=imag(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
        A_Matrix(j+3*N,2*k+6*N)=real(-Cji(k,j)*(-1)^j*(R2/(x1-x2))^k*(R1/(x1-x2))^j);
    end
    % (4N+1) 5N hang
    A_Matrix(j+4*N,2*j-1+2*N)=1;  
%     for k=1:N
%         A_Matrix(j+4*N,2*k-1+4*N)=-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j;
%     end
    for k=1:N
        A_Matrix(j+4*N,2*k-1+4*N)=real(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
        A_Matrix(j+4*N,2*k+4*N)=-imag(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
    end
    A_Matrix(j+4*N,2*j-1+6*N)=1;
    % (5N+1) 6N hang
    A_Matrix(j+5*N,2*j+2*N)=1;  
%     for k=1:N
%         A_Matrix(j+5*N,2*k+4*N)=-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j;
%     end
    for k=1:N
        A_Matrix(j+5*N,2*k-1+4*N)=imag(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
        A_Matrix(j+5*N,2*k+4*N)=real(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
    end
    A_Matrix(j+5*N,2*j+6*N)=-1;  
    % (6N+1) 7N hang
    A_Matrix(j+6*N,2*j-1+2*N)=gama20;  
%     for k=1:N
%         A_Matrix(j+6*N,2*k-1+4*N)=-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j;
%     end
    for k=1:N
        A_Matrix(j+6*N,2*k-1+4*N)=real(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
        A_Matrix(j+6*N,2*k+4*N)=-imag(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
    end
    A_Matrix(j+6*N,2*j-1+6*N)=-1;
    % (7N+1) 8N hang
    A_Matrix(j+7*N,2*j+2*N)=gama20;  
%     for k=1:N
%         A_Matrix(j+7*N,2*k+4*N)=-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j;
%     end
    for k=1:N
        A_Matrix(j+7*N,2*k-1+4*N)=imag(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
        A_Matrix(j+7*N,2*k+4*N)=real(-Cji(k,j)*(-1)^j*(R1/(x2-x1))^k*(R2/(x2-x1))^j);
    end
    A_Matrix(j+7*N,2*j+6*N)=1;     
end

end