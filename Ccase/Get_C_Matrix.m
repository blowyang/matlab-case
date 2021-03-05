function C_Matrix=Get_C_Matrix(K,N,a,R1)
C_Matrix=zeros(4*K,2*N);
for k=1:K
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    for n=1:N
        C_Matrix(k,2*n-1)=imag((R1/zk_AB)^n)-imag((R1/zk_CD)^n);
        C_Matrix(k,2*n)=real((R1/zk_AB)^n)-real((R1/zk_CD)^n);

        C_Matrix(k+K,2*n-1)=imag((R1/zk_DB)^n)-imag((R1/zk_CA)^n);
        C_Matrix(k+K,2*n)=real((R1/zk_DB)^n)-real((R1/zk_CA)^n);
        
        C_Matrix(k+2*K,2*n-1)=-real(n/R1*(R1/zk_AB)^(n+1))+real(n/R1*(R1/zk_CD)^(n+1));
        C_Matrix(k+2*K,2*n)=imag(n/R1*(R1/zk_AB)^(n+1))-imag(n/R1*(R1/zk_CD)^(n+1));
        
        C_Matrix(k+3*K,2*n-1)=-imag(n/R1*(R1/zk_DB)^(n+1))+imag(n/R1*(R1/zk_CA)^(n+1));
        C_Matrix(k+3*K,2*n)=-real(n/R1*(R1/zk_DB)^(n+1))+real(n/R1*(R1/zk_CA)^(n+1));
    end   
end
end