function [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm)
C_Matrix=zeros(8*K,4*N);
D_Matrix=zeros(8*K,4*M);
S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end
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

    [L1_zk_AB,L2_zk_AB,L3_zk_AB,L4_zk_AB,L5_zk_AB,L6_zk_AB]=Get_L_List(N,M,a,R,zk_AB);
    [L1_zk_CD,L2_zk_CD,L3_zk_CD,L4_zk_CD,L5_zk_CD,L6_zk_CD]=Get_L_List(N,M,a,R,zk_CD);
    [L1_zk_DB,L2_zk_DB,L3_zk_DB,L4_zk_DB,L5_zk_DB,L6_zk_DB]=Get_L_List(N,M,a,R,zk_DB);
    [L1_zk_CA,L2_zk_CA,L3_zk_CA,L4_zk_CA,L5_zk_CA,L6_zk_CA]=Get_L_List(N,M,a,R,zk_CA);
    %C_Matrix
    t_AB=[(km*L1_zk_AB.'*S1-zk_AB*conj(L3_zk_AB.'*S1)),-conj(L1_zk_AB.'*S1)];
    t_CD=[(km*L1_zk_CD.'*S1-zk_CD*conj(L3_zk_CD.'*S1)),-conj(L1_zk_CD.'*S1)];
    C_Matrix(k,:)=real(t_AB-t_CD)/2/Gm;
    C_Matrix(k+K,:)=imag(t_AB-t_CD)/2/Gm;  
    
    t_DB=[(km*L1_zk_DB.'*S1-zk_DB*conj(L3_zk_DB.'*S1)),-conj(L1_zk_DB.'*S1)];
    t_CA=[(km*L1_zk_CA.'*S1-zk_CA*conj(L3_zk_CA.'*S1)),-conj(L1_zk_CA.'*S1)];
    C_Matrix(k+2*K,:)=real(t_DB-t_CA)/2/Gm;
    C_Matrix(k+3*K,:)=imag(t_DB-t_CA)/2/Gm;
    
    t_AB=[(km*L3_zk_AB.'*S1-conj(L3_zk_AB.'*S1)+zk_AB*conj(L5_zk_AB.'*S1)),conj(L3_zk_AB.'*S1)];
    t_CD=[(km*L3_zk_CD.'*S1-conj(L3_zk_CD.'*S1)+zk_CD*conj(L5_zk_CD.'*S1)),conj(L3_zk_CD.'*S1)];
    C_Matrix(k+4*K,:)=real(t_AB-t_CD);
    C_Matrix(k+5*K,:)=imag(t_AB-t_CD);
    
    t_DB=[(km*L3_zk_DB.'*S1-conj(L3_zk_DB.'*S1)-zk_DB*conj(L5_zk_DB.'*S1)),-conj(L3_zk_DB.'*S1)];
    t_CA=[(km*L3_zk_CA.'*S1-conj(L3_zk_CA.'*S1)-zk_CA*conj(L5_zk_CA.'*S1)),-conj(L3_zk_CA.'*S1)];
    C_Matrix(k+6*K,:)=real(t_DB-t_CA);
    C_Matrix(k+7*K,:)=imag(t_DB-t_CA);
    %D_Matrix
    t_AB=[(km*L2_zk_AB.'*S0-zk_AB*conj(L4_zk_AB.'*S0)),-conj(L2_zk_AB.'*S0)];
    t_CD=[(km*L2_zk_CD.'*S0-zk_CD*conj(L4_zk_CD.'*S0)),-conj(L2_zk_CD.'*S0)];
    D_Matrix(k,:)=real(t_AB-t_CD)/2/Gm;
    D_Matrix(k+K,:)=imag(t_AB-t_CD)/2/Gm;  
    
    t_DB=[(km*L2_zk_DB.'*S0-zk_DB*conj(L4_zk_DB.'*S0)),-conj(L2_zk_DB.'*S0)];
    t_CA=[(km*L2_zk_CA.'*S0-zk_CA*conj(L4_zk_CA.'*S0)),-conj(L2_zk_CA.'*S0)];
    D_Matrix(k+2*K,:)=real(t_DB-t_CA)/2/Gm;
    D_Matrix(k+3*K,:)=imag(t_DB-t_CA)/2/Gm;
    
    t_AB=[(km*L4_zk_AB.'*S0-conj(L4_zk_AB.'*S0)+zk_AB*conj(L6_zk_AB.'*S0)),conj(L4_zk_AB.'*S0)];
    t_CD=[(km*L4_zk_CD.'*S0-conj(L4_zk_CD.'*S0)+zk_CD*conj(L6_zk_CD.'*S0)),conj(L4_zk_CD.'*S0)];
    D_Matrix(k+4*K,:)=real(t_AB-t_CD);
    D_Matrix(k+5*K,:)=imag(t_AB-t_CD);
    
    t_DB=[(km*L4_zk_DB.'*S0-conj(L4_zk_DB.'*S0)-zk_DB*conj(L6_zk_DB.'*S0)),-conj(L4_zk_DB.'*S0)];
    t_CA=[(km*L4_zk_CA.'*S0-conj(L4_zk_CA.'*S0)-zk_CA*conj(L6_zk_CA.'*S0)),-conj(L4_zk_CA.'*S0)];
    D_Matrix(k+6*K,:)=real(t_DB-t_CA);
    D_Matrix(k+7*K,:)=imag(t_DB-t_CA);
end
end