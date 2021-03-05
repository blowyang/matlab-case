function D_Matrix=Get_D_Matrix(K,M,a)
D_Matrix=zeros(4*K,2*M);
for k=1:K
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    for m=1:M
        D_Matrix(k,2*m-1)=imag(Get_P(m,a,zk_AB))-imag(Get_P(m,a,zk_CD));
        D_Matrix(k,2*m)=real(Get_P(m,a,zk_AB))-real(Get_P(m,a,zk_CD));

        D_Matrix(k+K,2*m-1)=imag(Get_P(m,a,zk_DB))-imag(Get_P(m,a,zk_CA));
        D_Matrix(k+K,2*m)=real(Get_P(m,a,zk_DB))-real(Get_P(m,a,zk_CA));
        
        D_Matrix(k+2*K,2*m-1)=real(Get_P_Diff(m,a,zk_AB))-real(Get_P_Diff(m,a,zk_CD));
        D_Matrix(k+2*K,2*m)=-imag(Get_P_Diff(m,a,zk_AB))+imag(Get_P_Diff(m,a,zk_CD));
        
        D_Matrix(k+3*K,2*m-1)=imag(Get_P_Diff(m,a,zk_DB))-imag(Get_P_Diff(m,a,zk_CA));
        D_Matrix(k+3*K,2*m)=real(Get_P_Diff(m,a,zk_DB))-real(Get_P_Diff(m,a,zk_CA));
    end   
end
end