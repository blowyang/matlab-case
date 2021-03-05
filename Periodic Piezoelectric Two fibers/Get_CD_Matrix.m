function [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,R1,R2,x1,x2,a,km,Gm)
C_Matrix=zeros(8*K,8*N);
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

    [L11_zk_AB,L12_zk_AB,L2_zk_AB,L31_zk_AB,L32_zk_AB,...
        L4_zk_AB,L51_zk_AB,L52_zk_AB,L6_zk_AB,~,~,~,~]=Get_L_List(N,M,R1,R2,x1,x2,a,zk_AB);

    [L11_zk_CD,L12_zk_CD,L2_zk_CD,L31_zk_CD,L32_zk_CD,...
        L4_zk_CD,L51_zk_CD,L52_zk_CD,L6_zk_CD,~,~,~,~]=Get_L_List(N,M,R1,R2,x1,x2,a,zk_CD);

    [L11_zk_DB,L12_zk_DB,L2_zk_DB,L31_zk_DB,L32_zk_DB,...
        L4_zk_DB,L51_zk_DB,L52_zk_DB,L6_zk_DB,~,~,~,~]=Get_L_List(N,M,R1,R2,x1,x2,a,zk_DB);

    [L11_zk_CA,L12_zk_CA,L2_zk_CA,L31_zk_CA,L32_zk_CA,...
        L4_zk_CA,L51_zk_CA,L52_zk_CA,L6_zk_CA,~,~,~,~]=Get_L_List(N,M,R1,R2,x1,x2,a,zk_CA);
    %C_Matrix
    t_AB=[(km*L11_zk_AB.'*S1-zk_AB*conj(L31_zk_AB.'*S1)),(km*L12_zk_AB.'*S1-zk_AB*conj(L32_zk_AB.'*S1)),...
        -conj(L11_zk_AB.'*S1),-conj(L12_zk_AB.'*S1)];
    t_CD=[(km*L11_zk_CD.'*S1-zk_CD*conj(L31_zk_CD.'*S1)),(km*L12_zk_CD.'*S1-zk_CD*conj(L32_zk_CD.'*S1)),...
        -conj(L11_zk_CD.'*S1),-conj(L12_zk_CD.'*S1)];
    C_Matrix(k,:)=real(t_AB-t_CD)/2/Gm;
    C_Matrix(k+K,:)=imag(t_AB-t_CD)/2/Gm;  
    
    t_DB=[(km*L11_zk_DB.'*S1-zk_DB*conj(L31_zk_DB.'*S1)),(km*L12_zk_DB.'*S1-zk_DB*conj(L32_zk_DB.'*S1)),...
        -conj(L11_zk_DB.'*S1),-conj(L12_zk_DB.'*S1)];
    t_CA=[(km*L11_zk_CA.'*S1-zk_CA*conj(L31_zk_CA.'*S1)),(km*L12_zk_CA.'*S1-zk_CA*conj(L32_zk_CA.'*S1)),...
        -conj(L11_zk_CA.'*S1),-conj(L12_zk_CA.'*S1)];
    C_Matrix(k+2*K,:)=real(t_DB-t_CA)/2/Gm;
    C_Matrix(k+3*K,:)=imag(t_DB-t_CA)/2/Gm;  
    
    t_AB=[km*L31_zk_AB.'*S1-conj(L31_zk_AB.'*S1)+zk_AB*conj(L51_zk_AB.'*S1),...
        km*L32_zk_AB.'*S1-conj(L32_zk_AB.'*S1)+zk_AB*conj(L52_zk_AB.'*S1),...
        conj(L31_zk_AB.'*S1),conj(L32_zk_AB.'*S1)];
    t_CD=[km*L31_zk_CD.'*S1-conj(L31_zk_CD.'*S1)+zk_CD*conj(L51_zk_CD.'*S1),...
        km*L32_zk_CD.'*S1-conj(L32_zk_CD.'*S1)+zk_CD*conj(L52_zk_CD.'*S1),...
        conj(L31_zk_CD.'*S1),conj(L32_zk_CD.'*S1)];
    C_Matrix(k+4*K,:)=real(t_AB-t_CD);
    C_Matrix(k+5*K,:)=imag(t_AB-t_CD);  

    
    
    t_DB=[km*L31_zk_DB.'*S1-conj(L31_zk_DB.'*S1)-zk_DB*conj(L51_zk_DB.'*S1),...
        km*L32_zk_DB.'*S1-conj(L32_zk_DB.'*S1)-zk_DB*conj(L52_zk_DB.'*S1),...
        -conj(L31_zk_DB.'*S1),-conj(L32_zk_DB.'*S1)];
    t_CA=[km*L31_zk_CA.'*S1-conj(L31_zk_CA.'*S1)-zk_CA*conj(L51_zk_CA.'*S1),...
        km*L32_zk_CA.'*S1-conj(L32_zk_CA.'*S1)-zk_CA*conj(L52_zk_CA.'*S1),...
        -conj(L31_zk_CA.'*S1),-conj(L32_zk_CA.'*S1)];
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
    
    t_AB=[km*L4_zk_AB.'*S0-conj(L4_zk_AB.'*S0)+zk_AB*conj(L6_zk_AB.'*S0),conj(L4_zk_AB.'*S0)];
    t_CD=[km*L4_zk_CD.'*S0-conj(L4_zk_CD.'*S0)+zk_CD*conj(L6_zk_CD.'*S0),conj(L4_zk_CD.'*S0)];
    D_Matrix(k+4*K,:)=real(t_AB-t_CD);
    D_Matrix(k+5*K,:)=imag(t_AB-t_CD);
    
    t_DB=[km*L4_zk_DB.'*S0-conj(L4_zk_DB.'*S0)-zk_DB*conj(L6_zk_DB.'*S0),-conj(L4_zk_DB.'*S0)];
    t_CA=[km*L4_zk_CA.'*S0-conj(L4_zk_CA.'*S0)-zk_CA*conj(L6_zk_CA.'*S0),-conj(L4_zk_CA.'*S0)];
    D_Matrix(k+6*K,:)=real(t_DB-t_CA);
    D_Matrix(k+7*K,:)=imag(t_DB-t_CA);
end
end