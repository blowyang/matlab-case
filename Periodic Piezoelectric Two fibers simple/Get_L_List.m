function [L11,L12,L2,L31,L32,L4,L51,L52,L6,L71,L72,L81,L82]=Get_L_List(N,R1,R2,x1,x2,a,z)
L11=zeros(N,1);
L12=zeros(N,1);
L31=zeros(N,1);
L32=zeros(N,1);
L51=zeros(N,1);
L52=zeros(N,1);
L71=zeros(N,1);
L72=zeros(N,1);
L81=zeros(N,1);
L82=zeros(N,1);

L2=zeros(N,1);
L4=zeros(N,1);
L6=zeros(N,1);

for n=1:N
    L11(n,1)=(R1/(z-x1))^n;
    L12(n,1)=(R2/(z-x2))^n;
    L31(n,1)=-n/R1*(R1/(z-x1))^(n+1);
    L32(n,1)=-n/R2*(R2/(z-x2))^(n+1);
    L51(n,1)=n*(n+1)/R1^2*(R1/(z-x1))^(n+2);
    L52(n,1)=n*(n+1)/R2^2*(R2/(z-x2))^(n+2);
    L71(n,1)=((z-x1)/R1)^n;
    L72(n,1)=((z-x2)/R2)^n;
    L81(n,1)=n/R1*((z-x1)/R1)^(n-1);
    L82(n,1)=n/R2*((z-x2)/R2)^(n-1);
end
for m=1:N
    L2(m,1)=(z/a)^m;
    L4(m,1)=m/a*(z/a)^(m-1);
    L6(m,1)=m*(m-1)/a^2*(z/a)^(m-2);
end
end