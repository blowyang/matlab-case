function C_Matrix=Get_C_Matrix(K,N,a,R1,R2,x1,x2)
C_Matrix=zeros(4*K,4*N);
S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end

for k=1:K
    zk_AB=k*a/(K+1)-a*(1-1i)/2;
    zk_CD=conj(zk_AB);
    zk_DB=k*a*1i/(K+1)+a*(1-1i)/2;
    zk_CA=conj(-zk_DB);
    
    L1_DB=zeros(N,1);
    L2_DB=zeros(N,1);
    L3_DB=zeros(N,1);
    L4_DB=zeros(N,1);
    L5_DB=zeros(N,1);
    L6_DB=zeros(N,1);
    
    L1_CA=zeros(N,1);
    L2_CA=zeros(N,1);
    L3_CA=zeros(N,1);
    L4_CA=zeros(N,1);
    L5_CA=zeros(N,1);
    L6_CA=zeros(N,1);
    
    L1_AB=zeros(N,1);
    L2_AB=zeros(N,1);
    L3_AB=zeros(N,1);
    L4_AB=zeros(N,1);
    L5_AB=zeros(N,1);
    L6_AB=zeros(N,1);
    
    L1_CD=zeros(N,1);
    L2_CD=zeros(N,1);
    L3_CD=zeros(N,1);
    L4_CD=zeros(N,1);
    L5_CD=zeros(N,1);
    L6_CD=zeros(N,1);
    
    for n=1:N
        L1_DB(n,1)=(R1/(zk_DB-x1))^n;
        L2_DB(n,1)=(R2/(zk_DB-x2))^n;
        L3_DB(n,1)=Get_P(n,a,zk_DB);
        L4_DB(n,1)=-n/R1*(R1/(zk_DB-x1))^(n+1);
        L5_DB(n,1)=-n/R2*(R2/(zk_DB-x2))^(n+1);
        L6_DB(n,1)=Get_P_Diff(n,a,zk_DB);
        
        L1_CA(n,1)=(R1/(zk_CA-x1))^n;
        L2_CA(n,1)=(R2/(zk_CA-x2))^n;
        L3_CA(n,1)=Get_P(n,a,zk_CA);
        L4_CA(n,1)=-n/R1*(R1/(zk_CA-x1))^(n+1);
        L5_CA(n,1)=-n/R2*(R2/(zk_CA-x2))^(n+1);
        L6_CA(n,1)=Get_P_Diff(n,a,zk_CA);
        
        L1_AB(n,1)=(R1/(zk_AB-x1))^n;
        L2_AB(n,1)=(R2/(zk_AB-x2))^n;
        L3_AB(n,1)=Get_P(n,a,zk_AB);
        L4_AB(n,1)=-n/R1*(R1/(zk_AB-x1))^(n+1);
        L5_AB(n,1)=-n/R2*(R2/(zk_AB-x2))^(n+1);
        L6_AB(n,1)=Get_P_Diff(n,a,zk_AB);
        
        L1_CD(n,1)=(R1/(zk_CD-x1))^n;
        L2_CD(n,1)=(R2/(zk_CD-x2))^n;
        L3_CD(n,1)=Get_P(n,a,zk_CD);
        L4_CD(n,1)=-n/R1*(R1/(zk_CD-x1))^(n+1);
        L5_CD(n,1)=-n/R2*(R2/(zk_CD-x2))^(n+1);
        L6_CD(n,1)=Get_P_Diff(n,a,zk_CD);
    end
    t_DB=[L1_DB.'*S1,L2_DB.'*S1];
    t_CA=[L1_CA.'*S1,L2_CA.'*S1];
    C_Matrix(k,:)=imag(t_DB-t_CA);
    
    t_AB=[L1_AB.'*S1,L2_AB.'*S1];
    t_CD=[L1_CD.'*S1,L2_CD.'*S1];
    C_Matrix(k+K,:)=imag(t_AB-t_CD);
    
    t_DB=[L4_DB.'*S1,L5_DB.'*S1];
    t_CA=[L4_CA.'*S1,L5_CA.'*S1];
    C_Matrix(k+2*K,:)=imag(t_DB-t_CA);
    
    t_AB=[L4_AB.'*S1,L5_AB.'*S1];
    t_CD=[L4_CD.'*S1,L5_CD.'*S1];
    C_Matrix(k+3*K,:)=real(t_AB-t_CD);
end
end