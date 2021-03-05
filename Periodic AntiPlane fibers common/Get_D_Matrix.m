function D_Matrix=Get_D_Matrix(K,M,a)
D_Matrix=zeros(4*K,2*M);

S0=zeros(M,2*M);
for m=1:M
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i;
end

for k=1:K
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    
    
    L3_DB=zeros(M,1);   
    L6_DB=zeros(M,1);
    
    L3_CA=zeros(M,1);
    L6_CA=zeros(M,1);

    L3_AB=zeros(M,1); 
    L6_AB=zeros(M,1);
    
    L3_CD=zeros(M,1); 
    L6_CD=zeros(M,1);
    
    for m=1:M
     
        L3_DB(m,1)=Get_P(m,a,zk_DB);
        L6_DB(m,1)=Get_P_Diff(m,a,zk_DB);

        L3_CA(m,1)=Get_P(m,a,zk_CA);
        L6_CA(m,1)=Get_P_Diff(m,a,zk_CA);
        
        L3_AB(m,1)=Get_P(m,a,zk_AB);
        L6_AB(m,1)=Get_P_Diff(m,a,zk_AB);
        
        L3_CD(m,1)=Get_P(m,a,zk_CD);
        L6_CD(m,1)=Get_P_Diff(m,a,zk_CD);
        
    end
    t_DB=L3_DB.'*S0;
    t_CA=L3_CA.'*S0;
    D_Matrix(k,:)=imag(t_DB-t_CA);
    
    t_AB=L3_AB.'*S0;
    t_CD=L3_CD.'*S0;
    D_Matrix(k+K,:)=imag(t_AB-t_CD);
    
    t_DB=L6_DB.'*S0;
    t_CA=L6_CA.'*S0;
    D_Matrix(k+2*K,:)=imag(t_DB-t_CA);
    
    t_AB=L6_AB.'*S0;
    t_CD=L6_CD.'*S0;
    D_Matrix(k+3*K,:)=real(t_AB-t_CD);
end
end