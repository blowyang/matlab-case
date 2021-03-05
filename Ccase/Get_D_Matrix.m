function D_Matrix=Get_D_Matrix(K,N,a)
D_Matrix=zeros(4*K,2*N);
for k=1:K
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    for m=1:N
        D_Matrix(k,2*m-1)=imag((zk_AB/a)^m)-imag((zk_CD/a)^m);
        D_Matrix(k,2*m)=real((zk_AB/a)^m)-real((zk_CD/a)^m);

        D_Matrix(k+K,2*m-1)=imag((zk_DB/a)^m)-imag((zk_CA/a)^m);
        D_Matrix(k+K,2*m)=real((zk_DB/a)^m)-real((zk_CA/a)^m);
        
        D_Matrix(k+2*K,2*m-1)=real(m/a*(zk_AB/a)^(m-1))-real(m/a*(zk_CD/a)^(m-1));
        D_Matrix(k+2*K,2*m)=-imag(m/a*(zk_AB/a)^(m-1))+imag(m/a*(zk_AB/a)^(m-1));
        
        D_Matrix(k+3*K,2*m-1)=imag(m/a*(zk_DB/a)^(m-1))-imag(m/a*(zk_CA/a)^(m-1));
        D_Matrix(k+3*K,2*m)=real(m/a*(zk_DB/a)^(m-1))-real(m/a*(zk_CA/a)^(m-1));
    end   
end
end