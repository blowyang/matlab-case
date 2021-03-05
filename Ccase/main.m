% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.

R2=1;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=1250;
gama1=0.1;
G1=G0/gama1;
gama2=0.1;
G2=G1/gama2;
R1=R2+0.1;
lambda=R2/R1;
lambda1=R1/a ;

N=6;
K=4;
[E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
B_Matrix=Get_B_Matrix(N,lambda1,gama1);

C_Matrix=Get_C_Matrix(K,N,a,R1);
D_Matrix=Get_D_Matrix(K,N,a);

H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

P1=zeros(4*K,1);
P2=zeros(4*K,1);

for i=1:K
    P1(i,1)=1;
    P2(i+K,1)=1;
end

%求delt
Xe1=((H_Matrix'*H_Matrix))\H_Matrix'*P1;
Xe2=((H_Matrix'*H_Matrix))\H_Matrix'*P2;

Xd1=-N_Matrix*B_Matrix*Xe1;
Xd2=-N_Matrix*B_Matrix*Xe2;

NP=100;
del=a/(NP+1);

S1=zeros(N,2*N);
S2=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end
for m=1:N
    S2(m,2*m-1)=1;
    S2(m,2*m)=1i; 
end
integral1=0;
integral2=0;
x_cof=zeros(2,2);
for np=1:NP
    z_AB=del*np-a/2+a*1i/2;
    L1=zeros(N,1);
    L2=zeros(N,1);
    for n=1:N
        L1(n,1)=n*(R1/z_AB)^(n+1);
    end
    for m=1:N
        L2(m,1)=Get_P_Diff(m,a,z_AB);
    end    
    integral1=integral1+real(Xe2'*S2'*L2-1/R1*Xd2'*S1'*L1)*del;
    integral2=integral2+real(Xe1'*S2'*L2-1/R1*Xd1'*S1'*L1)*del;
end
x_cof(1,1)=integral1;
x_cof(1,2)=integral2;
integral1=0;
integral2=0;
for np=1:NP
    z_DB=del*np*1i-a*1i/2+a/2;
    L1=zeros(N,1);
    L2=zeros(N,1);
    for n=1:N
        L1(n,1)=n*(R1/z_DB)^(n+1);
    end
    for m=1:N
        L2(m,1)=Get_P_Diff(m,a,z_DB);
    end    
    integral1=integral1+imag(Xe2'*S2'*L2-1/R1*Xd2'*S1'*L1)*del;
    integral2=integral2+imag(Xe1'*S2'*L2-1/R1*Xd1'*S1'*L1)*del;
end
x_cof(2,1)=integral1;
x_cof(2,2)=integral2;
y_res(1,1)=0.1;
y_res(2,1)=0.1;
res=x_cof\y_res;
Xd=res(1,1)*a*Xd2+res(2,1)*a*Xd1;
Xe=res(1,1)*a*Xe2+res(2,1)*a*Xe1;
w_AB=zeros(NP,1);
w_CD=zeros(NP,1);
w_AB_Sub_CD=zeros(NP,1);
w_DB=zeros(NP,1);
w_CA=zeros(NP,1);
w_DB_Sub_CA=zeros(NP,1);
S23_AB=zeros(NP,1);
S13_DB=zeros(NP,1);
for np=1:NP
    z_AB=del*np-a/2+a*1i/2;
    z_DB=del*np*1i-a*1i/2+a/2;
    L1=zeros(N,1);
    L2=zeros(N,1);
    for n=1:N
        L1(n,1)=(R1/z_AB)^n;
    end
    for m=1:N
        L2(m,1)=Get_P(m,a,z_AB);
    end   
    w_AB(np,1)=imag(Xd'*S1'*L1+Xe'*S2'*L2);
    for n=1:N
        L1(n,1)=n*(R1/z_AB)^n;
    end
    for m=1:N
        L2(m,1)=Get_P_Diff(m,a,z_AB);
    end
    S23_AB(np,1)=real(-1/R1*Xd'*S1'*L1+Xe'*S2'*L2);
    z_CD=conj(z_AB);
    for n=1:N
        L1(n,1)=(R1/z_CD)^n;
    end
    for m=1:N
        L2(m,1)=Get_P(m,a,z_CD);
    end   
    w_CD(np,1)=imag(Xd'*S1'*L1+Xe'*S2'*L2);
    w_AB_Sub_CD(np,1)=(w_AB(np,1)-w_CD(np,1))/(res(2,1)*a);
    
    for n=1:N
        L1(n,1)=(R1/z_DB)^n;
    end
    for m=1:N
        L2(m,1)=Get_P(m,a,z_DB);
    end   
    w_DB(np,1)=imag(Xd'*S1'*L1+Xe'*S2'*L2);
    for n=1:N
        L1(n,1)=n*(R1/z_DB)^n;
    end
    for m=1:N
        L2(m,1)=Get_P_Diff(m,a,z_DB);
    end
    S13_DB(np,1)=imag(-1/R1*Xd'*S1'*L1+Xe'*S2'*L2);
    
    z_CA=conj(-z_DB);
    for n=1:N
        L1(n,1)=(R1/z_CA)^n;
    end
    for m=1:N
        L2(m,1)=Get_P(m,a,z_CA);
    end   
    w_CA(np,1)=imag(Xd'*S1'*L1+Xe'*S2'*L2);
    w_DB_Sub_CA(np,1)=(w_DB(np,1)-w_CA(np,1))/(res(1,1)*a);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
