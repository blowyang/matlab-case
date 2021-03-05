function [delta1,delta2,delta3,delta4,delta_cof]=Get_Delta(a,R,Xc0,Xc1,Xc2,Xc3,Xc4,Xd0,Xd1,Xd2,Xd3,Xd4,...
    Xe0,Xe1,Xe2,Xe3,Xe4,Xf0,Xf1,Xf2,Xf3,Xf4,S022,S011,S012)
%delta_cof is equal module /G0
N=size(Xc1,1)/2;

NP=100;
del=a/(NP+1);

S=zeros(N,2*N);
for n=1:N
    S(n,2*n-1)=1;
    S(n,2*n)=1i;
end

delta_cof=zeros(4,4);
%
inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc1+L4_zk_AB.'*S*Xd1+conj(L3_zk_AB.'*S*Xc1+L4_zk_AB.'*S*Xd1)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc1+L6_zk_AB.'*S*Xd1)+conj(L3_zk_AB.'*S*Xe1+L4_zk_AB.'*S*Xf1);
    inte=inte+real(temp)*del;   
end
delta_cof(1,1)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc2+L4_zk_AB.'*S*Xd2+conj(L3_zk_AB.'*S*Xc2+L4_zk_AB.'*S*Xd2)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc2+L6_zk_AB.'*S*Xd2)+conj(L3_zk_AB.'*S*Xe2+L4_zk_AB.'*S*Xf2);
    inte=inte+real(temp)*del;   
end
delta_cof(1,2)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc3+L4_zk_AB.'*S*Xd3+conj(L3_zk_AB.'*S*Xc3+L4_zk_AB.'*S*Xd3)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc3+L6_zk_AB.'*S*Xd3)+conj(L3_zk_AB.'*S*Xe3+L4_zk_AB.'*S*Xf3);
    inte=inte+real(temp)*del;   
end
delta_cof(1,3)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc4+L4_zk_AB.'*S*Xd4+conj(L3_zk_AB.'*S*Xc4+L4_zk_AB.'*S*Xd4)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc4+L6_zk_AB.'*S*Xd4)+conj(L3_zk_AB.'*S*Xe4+L4_zk_AB.'*S*Xf4);
    inte=inte+real(temp)*del;   
end
delta_cof(1,4)=inte/a;
%
%
inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc1+L4_zk_AB.'*S*Xd1+conj(L3_zk_AB.'*S*Xc1+L4_zk_AB.'*S*Xd1)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc1+L6_zk_AB.'*S*Xd1)+conj(L3_zk_AB.'*S*Xe1+L4_zk_AB.'*S*Xf1);
    inte=inte+imag(temp)*del;   
end
delta_cof(2,1)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc2+L4_zk_AB.'*S*Xd2+conj(L3_zk_AB.'*S*Xc2+L4_zk_AB.'*S*Xd2)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc2+L6_zk_AB.'*S*Xd2)+conj(L3_zk_AB.'*S*Xe2+L4_zk_AB.'*S*Xf2);
    inte=inte+imag(temp)*del;   
end
delta_cof(2,2)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc3+L4_zk_AB.'*S*Xd3+conj(L3_zk_AB.'*S*Xc3+L4_zk_AB.'*S*Xd3)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc3+L6_zk_AB.'*S*Xd3)+conj(L3_zk_AB.'*S*Xe3+L4_zk_AB.'*S*Xf3);
    inte=inte+imag(temp)*del;   
end
delta_cof(2,3)=inte/a;

inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc4+L4_zk_AB.'*S*Xd4+conj(L3_zk_AB.'*S*Xc4+L4_zk_AB.'*S*Xd4)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc4+L6_zk_AB.'*S*Xd4)+conj(L3_zk_AB.'*S*Xe4+L4_zk_AB.'*S*Xf4);
    inte=inte+imag(temp)*del;   
end
delta_cof(2,4)=inte/a;
%
%
inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc1+L4_zk_DB.'*S*Xd1+conj(L3_zk_DB.'*S*Xc1+L4_zk_DB.'*S*Xd1)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc1+L6_zk_DB.'*S*Xd1)-conj(L3_zk_DB.'*S*Xe1+L4_zk_DB.'*S*Xf1);
    inte=inte+real(temp)*del;   
end
delta_cof(3,1)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc2+L4_zk_DB.'*S*Xd2+conj(L3_zk_DB.'*S*Xc2+L4_zk_DB.'*S*Xd2)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc2+L6_zk_DB.'*S*Xd2)-conj(L3_zk_DB.'*S*Xe2+L4_zk_DB.'*S*Xf2);
    inte=inte+real(temp)*del;   
end
delta_cof(3,2)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc3+L4_zk_DB.'*S*Xd3+conj(L3_zk_DB.'*S*Xc3+L4_zk_DB.'*S*Xd3)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc3+L6_zk_DB.'*S*Xd3)-conj(L3_zk_DB.'*S*Xe3+L4_zk_DB.'*S*Xf3);
    inte=inte+real(temp)*del;    
end
delta_cof(3,3)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc4+L4_zk_DB.'*S*Xd4+conj(L3_zk_DB.'*S*Xc4+L4_zk_DB.'*S*Xd4)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc4+L6_zk_DB.'*S*Xd4)-conj(L3_zk_DB.'*S*Xe4+L4_zk_DB.'*S*Xf4);
    inte=inte+real(temp)*del;    
end
delta_cof(3,4)=inte/a;
%
%
inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc1+L4_zk_DB.'*S*Xd1+conj(L3_zk_DB.'*S*Xc1+L4_zk_DB.'*S*Xd1)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc1+L6_zk_DB.'*S*Xd1)-conj(L3_zk_DB.'*S*Xe1+L4_zk_DB.'*S*Xf1);
    inte=inte+imag(-temp)*del;   
end
delta_cof(4,1)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc2+L4_zk_DB.'*S*Xd2+conj(L3_zk_DB.'*S*Xc2+L4_zk_DB.'*S*Xd2)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc2+L6_zk_DB.'*S*Xd2)-conj(L3_zk_DB.'*S*Xe2+L4_zk_DB.'*S*Xf2);
    inte=inte+imag(-temp)*del;   
end
delta_cof(4,2)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc3+L4_zk_DB.'*S*Xd3+conj(L3_zk_DB.'*S*Xc3+L4_zk_DB.'*S*Xd3)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc3+L6_zk_DB.'*S*Xd3)-conj(L3_zk_DB.'*S*Xe3+L4_zk_DB.'*S*Xf3);
    inte=inte+imag(-temp)*del;    
end
delta_cof(4,3)=inte/a;

inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc4+L4_zk_DB.'*S*Xd4+conj(L3_zk_DB.'*S*Xc4+L4_zk_DB.'*S*Xd4)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc4+L6_zk_DB.'*S*Xd4)-conj(L3_zk_DB.'*S*Xe4+L4_zk_DB.'*S*Xf4);
    inte=inte+imag(-temp)*del;    
end
delta_cof(4,4)=inte/a;
%
inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc0+L4_zk_AB.'*S*Xd0+conj(L3_zk_AB.'*S*Xc0+L4_zk_AB.'*S*Xd0)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc0+L6_zk_AB.'*S*Xd0)+conj(L3_zk_AB.'*S*Xe0+L4_zk_AB.'*S*Xf0);
    inte=inte+real(temp)*del;   
end
Stress_res(1,1)=S022-inte;
%
%
inte=0;
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    L1_zk_AB=zeros(N,1);
    L2_zk_AB=zeros(N,1);
    L3_zk_AB=zeros(N,1);
    L4_zk_AB=zeros(N,1);
    L5_zk_AB=zeros(N,1);
    L6_zk_AB=zeros(N,1);
    for n=1:N
        L1_zk_AB(n,1)=(R/zk_AB)^n;
        L2_zk_AB(n,1)=(zk_AB/a)^n;
        L3_zk_AB(n,1)=-n/R*(R/zk_AB)^(n+1);
        L4_zk_AB(n,1)=n/a*(zk_AB/a)^(n-1);
        L5_zk_AB(n,1)=n*(n+1)/R^2*(R/zk_AB)^(n+2);
        L6_zk_AB(n,1)=n*(n-1)/a^2*(zk_AB/a)^(n-2);
    end
    temp=L3_zk_AB.'*S*Xc0+L4_zk_AB.'*S*Xd0+conj(L3_zk_AB.'*S*Xc0+L4_zk_AB.'*S*Xd0)+...
        conj(zk_AB)*(L5_zk_AB.'*S*Xc0+L6_zk_AB.'*S*Xd0)+conj(L3_zk_AB.'*S*Xe0+L4_zk_AB.'*S*Xf0);
    inte=inte+imag(temp)*del;   
end
Stress_res(2,1)=S012-inte;
%
inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc0+L4_zk_DB.'*S*Xd0+conj(L3_zk_DB.'*S*Xc0+L4_zk_DB.'*S*Xd0)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc0+L6_zk_DB.'*S*Xd0)-conj(L3_zk_DB.'*S*Xe0+L4_zk_DB.'*S*Xf0);
    inte=inte+real(temp)*del;    
end
Stress_res(3,1)=S011-inte;
%
inte=0;
for np=1:NP
    zk_DB=del*np*1i-a*1i/2+a/2;
    L1_zk_DB=zeros(N,1);
    L2_zk_DB=zeros(N,1);
    L3_zk_DB=zeros(N,1);
    L4_zk_DB=zeros(N,1);
    L5_zk_DB=zeros(N,1);
    L6_zk_DB=zeros(N,1);
    
    for n=1:N
        L1_zk_DB(n,1)=(R/zk_DB)^n;
        L2_zk_DB(n,1)=(zk_DB/a)^n;
        L3_zk_DB(n,1)=-n/R*(R/zk_DB)^(n+1);
        L4_zk_DB(n,1)=n/a*(zk_DB/a)^(n-1);
        L5_zk_DB(n,1)=n*(n+1)/R^2*(R/zk_DB)^(n+2);
        L6_zk_DB(n,1)=n*(n-1)/a^2*(zk_DB/a)^(n-2);
    end
    temp=L3_zk_DB.'*S*Xc0+L4_zk_DB.'*S*Xd0+conj(L3_zk_DB.'*S*Xc0+L4_zk_DB.'*S*Xd0)-...
        conj(zk_DB)*(L5_zk_DB.'*S*Xc0+L6_zk_DB.'*S*Xd0)-conj(L3_zk_DB.'*S*Xe0+L4_zk_DB.'*S*Xf0);
    inte=inte+imag(-temp)*del;    
end
Stress_res(4,1)=S012-inte;

%res=delta_cof\Stress_res;
res=(delta_cof'*delta_cof)\(delta_cof')*Stress_res;
delta1=res(1,1);
delta2=res(2,1);
delta3=res(3,1);
delta4=res(4,1);
end