function [L1,L2,L3,L4,L5,L6]=Get_L_List(N,M,a,R,z)
L1=zeros(N,1);
L2=zeros(M,1);
L3=zeros(N,1);
L4=zeros(M,1);
L5=zeros(N,1);
L6=zeros(M,1);

for n=1:N
    L1(n,1)=(R/z)^n;
    L3(n,1)=-n/R*(R/z)^(n+1);
    L5(n,1)=n*(n+1)/R^2*(R/z)^(n+2);
end
for m=1:M
    L2(m,1)=Get_P(m,a,z);
    L4(m,1)=Get_P_Diff(m,a,z);
    L6(m,1)=Get_P_2Diff(m,a,z);
end
end