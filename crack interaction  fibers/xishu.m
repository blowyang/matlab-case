function [A,B]=xishu(N,Z1,Z2,Z3,Zd,R1,R2,R3,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef)

B1=(S11+S22)/4;
C1=0;
B2=(S22-S11)/2;
C2=S12;

um=Em/2/(1+vm);
uf=Ef/2/(1+vf);

CE=um/uf;
km=3-4*vm;
kf=3-4*vf;

gama=um*(by-1i*bx)/pi/(1+km);
%gama=0;
A(1,1)=0;
B(1)=0;

%压电纤维1与基体应力连续条件、位移连续条件，一次多项式；
j=1;
i=1;
%系数ak
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数bk
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R2/(Z1-Z2))^(k+1)*R1/R2;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R2/(Z1-Z2))^(k+1)*R1/R2);
    j=j+1;
    A(i,j)=-R1*k/R2*(R2/conj(Z1-Z2))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R1*k/R2*(R2/conj(Z1-Z2))^(k+1));
    j=j+1;
end
%系数ik
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R3/(Z1-Z3))^(k+1)*R1/R3;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R3/(Z1-Z3))^(k+1)*R1/R3);
    j=j+1;
    A(i,j)=-R1*k/R3*(R3/conj(Z1-Z3))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R1*k/R3*(R3/conj(Z1-Z3))^(k+1));
    j=j+1;
end
%系数ck
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数jk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ek
for k=1:1:N
    if k==1
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE*1;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=-(B1+1i*C1)*R1-gama*(R1/(Z1-Zd))-R1*(B1-1i*C1)-conj(gama)*(R1/conj(Z1-Zd))-e31*E3*R1;
B(i+1)=4*km/(1+km)*(-(B1+1i*C1)*R1-gama*(R1/(Z1-Zd)))-4/(1+km)*(-R1*(B1-1i*C1)-conj(gama)*(R1/conj(Z1-Zd)));

%压电纤维1与基体应力连续条件、位移连续条件，二次多项式；
j=1;
i=i+2;
%系数ak
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-Z1/R1;
       A(i+1,j)=-4/(1+km)*(-Z1/R1);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数bk
for k=1:1:N
    A(i,j)=Cji(k,2)*(R2/(Z1-Z2))^(k+2)*(R1/R2)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R2/(Z1-Z2))^(k+2)*(R1/R2)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ik
for k=1:1:N
    A(i,j)=Cji(k,2)*(R3/(Z1-Z3))^(k+2)*(R1/R3)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R3/(Z1-Z3))^(k+2)*(R1/R3)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ck
for k=1:1:N
    if k==2
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数jk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ek
for k=1:1:N
    if k==2
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=gama/2*(R1/(Z1-Zd))^2;
B(i+1)=4*km/(1+km)*(gama/2*(R1/(Z1-Zd))^2);

%压电纤维1与基体应力连续条件、位移连续条件，n次多项式；
i=i+2;

for n=3:1:N
    j=1;
    %系数ak
    for k=1:1:N
        if k==n-2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n-2);
           A(i+1,j)=4/(1+km)*(n-2);
            j=j+1;
        elseif k==n-1
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-Z1/R1*(n-1);
           A(i+1,j)=4/(1+km)*Z1/R1*(n-1);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;            
        end
    end
    %系数bk
    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R2/(Z1-Z2))^(k+n)*(R1/R2)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R2/(Z1-Z2))^(k+n)*(R1/R2)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数ik
    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R3/(Z1-Z3))^(k+n)*(R1/R3)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R3/(Z1-Z3))^(k+n)*(R1/R3)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数ck
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=1;
           A(i+1,j)=-4/(1+km);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数dk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数jk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数ek
    for k=1:1:N
        if k==n
            A(i,j)=-1;
           A(i+1,j)=-4*kf/(1+kf)*CE;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
     %系数fk
   for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end 
      %系数gk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
     end
       %系数hk
   for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
       %系数kk
   for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
        %系数lk
   for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
   end
    
   B(i)=gama*(-1)^n/n*(R1/(Z1-Zd))^n;
    B(i+1)=4*km/(1+km)*(gama*(-1)^n/n*(R1/(Z1-Zd))^n);
    i=i+2;
end

%压电纤维1与基体应力连续条件、位移连续条件，负一次多项式；
j=1;
%系数ak
for k=1:1:N
    if k==1
        A(i,j)=1;
       A(i+1,j)=4*km/(1+km);
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数bk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R1*k/R2*Cji(k+1,2)*(R2/conj(Z1-Z2))^(k+3)*(R1/R2)^2+Z1*k/R2*Cji(k+1,1)*(R2/conj(Z1-Z2))^(k+2)*R1/R2;%change
   A(i+1,j)=-4/(1+km)*(-R1*k/R2*Cji(k+1,2)*(R2/conj(Z1-Z2))^(k+3)*(R1/R2)^2+Z1*k/R2*Cji(k+1,1)*(R2/conj(Z1-Z2))^(k+2)*R1/R2);%
    j=j+1;
end
%系数ik
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R1*k/R3*Cji(k+1,2)*(R3/conj(Z1-Z3))^(k+3)*(R1/R3)^2+Z1*k/R3*Cji(k+1,1)*(R3/conj(Z1-Z3))^(k+2)*R1/R3;%change
   A(i+1,j)=-4/(1+km)*(-R1*k/R3*Cji(k+1,2)*(R3/conj(Z1-Z3))^(k+3)*(R1/R3)^2+Z1*k/R3*Cji(k+1,1)*(R3/conj(Z1-Z3))^(k+2)*R1/R3);%
    j=j+1;
end
%系数ck
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R2/conj(Z1-Z2))^(k+1)*R1/R2;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R2/conj(Z1-Z2))^(k+1)*R1/R2);
    j=j+1;
end
%系数jk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R3/conj(Z1-Z3))^(k+1)*R1/R3;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R3/conj(Z1-Z3))^(k+1)*R1/R3);
    j=j+1;
end
%系数ek
for k=1:1:N
    if k==3
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-3;
       A(i+1,j)=4/(1+kf)*CE*3;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数fk

for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=-conj(gama)*(R1/conj(Z1-Zd))^3+Z1*conj(gama)/R1*(R1/conj(Z1-Zd))^2-(B2-1i*C2)*R1-gama*R1/conj(Z1-Zd)-conj(gama)*Zd/R1*(R1/conj(Z1-Zd))^2;
B(i+1)=-4/(1+km)*(-conj(gama)*(R1/conj(Z1-Zd))^3+Z1*conj(gama)/R1*(R1/conj(Z1-Zd))^2-(B2-1i*C2)*R1-gama*R1/conj(Z1-Zd)-conj(gama)*Zd/R1*(R1/conj(Z1-Zd))^2);
    
%压电纤维1与基体应力连续条件、位移连续条件，负n次多项式；
i=i+2;

for n=2:1:N
    j=1;
    %系数ak
    for k=1:1:N
        if k==n
            A(i,j)=1;
           A(i+1,j)=4*km/(1+km);
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;            
        end
    end
    %系数bk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-(-1)^(n+1)*R1*k/R2*Cji(k+1,n+1)*(R2/conj(Z1-Z2))^(k+n+2)*(R1/R2)^(n+1)-(-1)^n*Z1*k/R2*Cji(k+1,n)*(R2/conj(Z1-Z2))^(k+n+1)*(R1/R2)^n;
       A(i+1,j)=-4/(1+km)*(-(-1)^(n+1)*R1*k/R2*Cji(k+1,n+1)*(R2/conj(Z1-Z2))^(k+n+2)*(R1/R2)^(n+1)-(-1)^n*Z1*k/R2*Cji(k+1,n)*(R2/conj(Z1-Z2))^(k+n+1)*(R1/R2)^n);
        j=j+1;
    end
    %系数ik
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-(-1)^(n+1)*R1*k/R3*Cji(k+1,n+1)*(R3/conj(Z1-Z3))^(k+n+2)*(R1/R3)^(n+1)-(-1)^n*Z1*k/R3*Cji(k+1,n)*(R3/conj(Z1-Z3))^(k+n+1)*(R1/R3)^n;
       A(i+1,j)=-4/(1+km)*(-(-1)^(n+1)*R1*k/R3*Cji(k+1,n+1)*(R3/conj(Z1-Z3))^(k+n+2)*(R1/R3)^(n+1)-(-1)^n*Z1*k/R3*Cji(k+1,n)*(R3/conj(Z1-Z3))^(k+n+1)*(R1/R3)^n);
        j=j+1;
    end 
    %系数ck
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数dk
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R2/conj(Z1-Z2))^(k+n)*(R1/R2)^n;
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R2/conj(Z1-Z2))^(k+n)*(R1/R2)^n);
        j=j+1;
    end
    %系数jk
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R3/conj(Z1-Z3))^(k+n)*(R1/R3)^n;
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R3/conj(Z1-Z3))^(k+n)*(R1/R3)^n);
        j=j+1;
    end
    %系数ek
    for k=1:1:N
        if k==n+2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n+2);
           A(i+1,j)=4/(1+kf)*CE*(n+2);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数fk
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-1;
           A(i+1,j)=4/(1+kf)*CE;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数gk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数hk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
     %系数kk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
     %系数lk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
      
    B(i)=(-1)^n*conj(gama)*(R1/conj(Z1-Zd))^(n+2)-(-1)^n*Z1*conj(gama)/R1*(R1/conj(Z1-Zd))^(n+1)+(-1)^n*gama/n*(R1/conj(Z1-Zd))^n+(-1)^n*conj(gama)*Zd/R1*(R1/conj(Z1-Zd))^(n+1);
    B(i+1)=-4/(1+km)*((-1)^n*conj(gama)*(R1/conj(Z1-Zd))^(n+2)-(-1)^n*Z1*conj(gama)/R1*(R1/conj(Z1-Zd))^(n+1)+(-1)^n*gama/n*(R1/conj(Z1-Zd))^n+(-1)^n*conj(gama)*Zd/R1*(R1/conj(Z1-Zd))^(n+1));
    i=i+2;
end
%压电纤维2与基体应力连续条件、位移连续条件；
%压电纤维2与基体应力连续条件、位移连续条件；
%压电纤维2与基体应力连续条件、位移连续条件；

j=1;

% %压电纤维2与基体应力连续条件、位移连续条件，一次多项式；
% j=1;
% i=i+4;
%系数ak
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R1/(Z2-Z1))^(k+1)*R2/R1;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R1/(Z2-Z1))^(k+1)*R2/R1);
    j=j+1;
    A(i,j)=-R2*k/R1*(R1/conj(Z2-Z1))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R2*k/R1*(R1/conj(Z2-Z1))^(k+1));
    j=j+1;
end
%系数bk

for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ik
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R3/(Z2-Z3))^(k+1)*R2/R3;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R3/(Z2-Z3))^(k+1)*R2/R3);
    j=j+1;
    A(i,j)=-R2*k/R3*(R3/conj(Z2-Z3))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R2*k/R3*(R3/conj(Z2-Z3))^(k+1));
    j=j+1;
end
%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数dk
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end    
end
%系数jk
for k=1:1:N
    
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;

end
%系数ek
for k=1:1:N
    
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;

end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    if k==1
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=-(B1+1i*C1)*R2-gama*(R2/(Z2-Zd))-R2*(B1-1i*C1)-conj(gama)*(R2/conj(Z2-Zd))-e31*E3*R2;
B(i+1)=4*km/(1+km)*(-(B1+1i*C1)*R2-gama*(R2/(Z2-Zd)))-4/(1+km)*(-R2*(B1-1i*C1)-conj(gama)*(R2/conj(Z2-Zd)));

%压电纤维2与基体应力连续条件、位移连续条件，二次多项式；
j=1;
i=i+2;
%系数ak
for k=1:1:N
    A(i,j)=Cji(k,2)*(R1/(Z2-Z1))^(k+2)*(R2/R1)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R1/(Z2-Z1))^(k+2)*(R2/R1)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数bk
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-Z2/R2;
       A(i+1,j)=4/(1+km)*Z2/R2;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数ik
for k=1:1:N
    A(i,j)=Cji(k,2)*(R3/(Z2-Z3))^(k+2)*(R2/R3)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R3/(Z2-Z3))^(k+2)*(R2/R3)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数dk
for k=1:1:N
    if k==2
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数jk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ek
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    if k==2
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=gama/2*(R2/(Z2-Zd))^2;
B(i+1)=4*km/(1+km)*(gama/2*(R2/(Z2-Zd))^2);

%压电纤维2与基体应力连续条件、位移连续条件，n次多项式；
i=i+2;
for n=3:1:N
    j=1;
    %系数ak

    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R1/(Z2-Z1))^(k+n)*(R2/R1)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R1/(Z2-Z1))^(k+n)*(R2/R1)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数bk

    for k=1:1:N
        if k==n-2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n-2);
           A(i+1,j)=4/(1+km)*(n-2);
            j=j+1;
        elseif k==n-1
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-Z2/R2*(n-1);
           A(i+1,j)=4/(1+km)*Z2/R2*(n-1);
            j=j+1; 
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;         
        end
    end
    %系数ik
    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R3/(Z2-Z3))^(k+n)*(R2/R3)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R3/(Z2-Z3))^(k+n)*(R2/R3)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
%系数ck
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数dk
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=1;
           A(i+1,j)=-4/(1+km);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数jk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数ek
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数fk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end 
    %系数gk
     for k=1:1:N
        if k==n
            A(i,j)=-1;
           A(i+1,j)=-4*kf/(1+kf)*CE;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
     end
     %系数hk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
         %系数kk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
     %系数lk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end

    B(i)=gama*(-1)^n/n*(R2/(Z2-Zd))^n;
    B(i+1)=4*km/(1+km)*(gama*(-1)^n/n*(R2/(Z2-Zd))^n);
    i=i+2;
end

%压电纤维2与基体应力连续条件、位移连续条件，负一次多项式；
j=1;
%系数ak
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R2*k/R1*Cji(k+1,2)*(R1/conj(Z2-Z1))^(k+3)*(R2/R1)^2+Z2*k/R1*Cji(k+1,1)*(R1/conj(Z2-Z1))^(k+2)*R2/R1;
   A(i+1,j)=-4/(1+km)*(-R2*k/R1*Cji(k+1,2)*(R1/conj(Z2-Z1))^(k+3)*(R2/R1)^2+Z2*k/R1*Cji(k+1,1)*(R1/conj(Z2-Z1))^(k+2)*R2/R1);
    j=j+1;

end
%系数bk
for k=1:1:N
    if k==1
        A(i,j)=1;
       A(i+1,j)=4*km/(1+km);
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数ik
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R2*k/R3*Cji(k+1,2)*(R3/conj(Z2-Z3))^(k+3)*(R2/R3)^2+Z2*k/R3*Cji(k+1,1)*(R3/conj(Z2-Z3))^(k+2)*R2/R3;
   A(i+1,j)=-4/(1+km)*(-R2*k/R3*Cji(k+1,2)*(R3/conj(Z2-Z3))^(k+3)*(R2/R3)^2+Z2*k/R3*Cji(k+1,1)*(R3/conj(Z2-Z3))^(k+2)*R2/R3);
    j=j+1;

end
%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R1/conj(Z2-Z1))^(k+1)*R2/R1;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R1/conj(Z2-Z1))^(k+1)*R2/R1);
    j=j+1;
end
%系数dk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数jk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R3/conj(Z2-Z3))^(k+1)*R2/R3;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R3/conj(Z2-Z3))^(k+1)*R2/R3);
    j=j+1;
end
%系数ek
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
  
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
%系数gk
for k=1:1:N
    if k==3
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-3;
       A(i+1,j)=4/(1+kf)*CE*3;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数hk
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE*1;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数kk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
B(i)=-conj(gama)*(R2/conj(Z2-Zd))^3+Z2*conj(gama)/R2*(R2/conj(Z2-Zd))^2-(B2-1i*C2)*R2-gama*R2/conj(Z2-Zd)-conj(gama)*Zd/R2*(R2/conj(Z2-Zd))^2;
B(i+1)=-4/(1+km)*(-conj(gama)*(R2/conj(Z2-Zd))^3+Z2*conj(gama)/R2*(R2/conj(Z2-Zd))^2-(B2-1i*C2)*R2-gama*R2/conj(Z2-Zd)-conj(gama)*Zd/R2*(R2/conj(Z2-Zd))^2);
%压电纤维2与基体应力连续条件、位移连续条件，负n次多项式；
i=i+2;

for n=2:1:N
    j=1;
    %系数ak
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-R2*k/R1*(-1)^(n+1)*Cji(k+1,n+1)*(R1/conj(Z2-Z1))^(k+n+2)*(R2/R1)^(n+1)-Z2*k/R1*(-1)^n*Cji(k+1,n)*(R1/conj(Z2-Z1))^(k+n+1)*(R2/R1)^n;
       A(i+1,j)=-4/(1+km)*(-R2*k/R1*(-1)^(n+1)*Cji(k+1,n+1)*(R1/conj(Z2-Z1))^(k+n+2)*(R2/R1)^(n+1)-Z2*k/R1*(-1)^n*Cji(k+1,n)*(R1/conj(Z2-Z1))^(k+n+1)*(R2/R1)^n);
        j=j+1;
    end
    %系数bk
    for k=1:1:N
        if k==n
            A(i,j)=1;
            A(i+1,j)=4*km/(1+km)*1;
            j=j+1;
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数ik
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-R2*k/R3*(-1)^(n+1)*Cji(k+1,n+1)*(R3/conj(Z2-Z3))^(k+n+2)*(R2/R3)^(n+1)-Z2*k/R3*(-1)^n*Cji(k+1,n)*(R3/conj(Z2-Z3))^(k+n+1)*(R2/R3)^n;
       A(i+1,j)=-4/(1+km)*(-R2*k/R3*(-1)^(n+1)*Cji(k+1,n+1)*(R3/conj(Z2-Z3))^(k+n+2)*(R2/R3)^(n+1)-Z2*k/R3*(-1)^n*Cji(k+1,n)*(R3/conj(Z2-Z3))^(k+n+1)*(R2/R3)^n);
        j=j+1;
    end
%系数ck
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R1/conj(Z2-Z1))^(k+n)*(R2/R1)^n;        
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R1/conj(Z2-Z1))^(k+n)*(R2/R1)^n);
        j=j+1;
    end
    %系数dk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数jk
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R3/conj(Z2-Z3))^(k+n)*(R2/R3)^n;        
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R3/conj(Z2-Z3))^(k+n)*(R2/R3)^n);
        j=j+1;
    end
    %系数ek
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数fk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end
    %系数gk
    for k=1:1:N
        if k==n+2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n+2);
           A(i+1,j)=4/(1+kf)*CE*(n+2);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数hk
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-1;
           A(i+1,j)=4/(1+kf)*CE*1;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数kk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end
    %系数lk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end
    
    B(i)=-conj(gama)*(-1)^(n+1)*(R2/conj(Z2-Zd))^(n+2)-Z2*conj(gama)*(-1)^n/R2*(R2/conj(Z2-Zd))^(n+1)+gama*(-1)^n/n*(R2/conj(Z2-Zd))^n+conj(gama)*(-1)^n*Zd/R2*(R2/conj(Z2-Zd))^(n+1);
    B(i+1)=-4/(1+km)*(-conj(gama)*(-1)^(n+1)*(R2/conj(Z2-Zd))^(n+2)-Z2*conj(gama)*(-1)^n/R2*(R2/conj(Z2-Zd))^(n+1)+gama*(-1)^n/n*(R2/conj(Z2-Zd))^n+conj(gama)*(-1)^n*Zd/R2*(R2/conj(Z2-Zd))^(n+1));
    i=i+2;
end
%压电纤维3与基体应力连续条件、位移连续条件；

j=1;

% %压电纤维3与基体应力连续条件、位移连续条件，一次多项式；
% j=1;
% i=i+4;
%系数ak
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R1/(Z3-Z1))^(k+1)*R3/R1;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R1/(Z3-Z1))^(k+1)*R3/R1);
    j=j+1;
    A(i,j)=-R3*k/R1*(R1/conj(Z3-Z1))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R3*k/R1*(R1/conj(Z3-Z1))^(k+1));
    j=j+1;
end
%系数bk
for k=1:1:N
    A(i,j)=-Cji(k,1)*(R2/(Z3-Z2))^(k+1)*R3/R2;
   A(i+1,j)=4*km/(1+km)*(-Cji(k,1)*(R2/(Z3-Z2))^(k+1)*R3/R2);
    j=j+1;
    A(i,j)=-R3*k/R2*(R2/conj(Z3-Z2))^(k+1);
   A(i+1,j)=-4/(1+km)*(-R3*k/R2*(R2/conj(Z3-Z2))^(k+1));
    j=j+1;
end
%系数ik
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end

%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数jk
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end    
end

%系数ek
for k=1:1:N
    
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;

end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end

%系数kk
for k=1:1:N
    if k==1
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
B(i)=-(B1+1i*C1)*R3-gama*(R3/(Z3-Zd))-R3*(B1-1i*C1)-conj(gama)*(R3/conj(Z3-Zd))-e31*E3*R3;
B(i+1)=4*km/(1+km)*(-(B1+1i*C1)*R3-gama*(R3/(Z3-Zd)))-4/(1+km)*(-R3*(B1-1i*C1)-conj(gama)*(R3/conj(Z3-Zd)));

%压电纤维3与基体应力连续条件、位移连续条件，二次多项式；
j=1;
i=i+2;
%系数ak
for k=1:1:N
    A(i,j)=Cji(k,2)*(R1/(Z3-Z1))^(k+2)*(R3/R1)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R1/(Z3-Z1))^(k+2)*(R3/R1)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数bk
for k=1:1:N
    A(i,j)=Cji(k,2)*(R2/(Z3-Z2))^(k+2)*(R3/R2)^2;
   A(i+1,j)=4*km/(1+km)*(Cji(k,2)*(R2/(Z3-Z2))^(k+2)*(R3/R2)^2);
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数ik
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-Z3/R3;
       A(i+1,j)=4/(1+km)*Z3/R3;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end

%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数jk
for k=1:1:N
    if k==2
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=1;
       A(i+1,j)=-4/(1+km);
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数ek
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end
%系数kk
for k=1:1:N
    if k==2
        A(i,j)=-1;
       A(i+1,j)=-4*kf/(1+kf)*CE;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end     
end
%系数lk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
end

B(i)=gama/2*(R3/(Z3-Zd))^2;
B(i+1)=4*km/(1+km)*(gama/2*(R3/(Z3-Zd))^2);

%压电纤维1与基体应力连续条件、位移连续条件，n次多项式；
i=i+2;
for n=3:1:N
    j=1;
    %系数ak

    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R1/(Z3-Z1))^(k+n)*(R3/R1)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R1/(Z3-Z1))^(k+n)*(R3/R1)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数bk
    for k=1:1:N
        A(i,j)=(-1)^n*Cji(k,n)*(R2/(Z3-Z2))^(k+n)*(R3/R2)^n;
       A(i+1,j)=4*km/(1+km)*(-1)^n*Cji(k,n)*(R2/(Z3-Z2))^(k+n)*(R3/R2)^n;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end

    %系数ik
    for k=1:1:N
        if k==n-2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n-2);
           A(i+1,j)=4/(1+km)*(n-2);
            j=j+1;
        elseif k==n-1
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-Z3/R3*(n-1);
           A(i+1,j)=4/(1+km)*Z3/R3*(n-1);
            j=j+1; 
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;         
        end
    end
%系数ck
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数dk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end

    %系数jk
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=1;
           A(i+1,j)=-4/(1+km);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数ek
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数fk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
         %系数gk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
     %系数hk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    
    %系数kk
     for k=1:1:N
        if k==n
            A(i,j)=-1;
           A(i+1,j)=-4*kf/(1+kf)*CE;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
     end
     %系数lk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end

    B(i)=gama*(-1)^n/n*(R3/(Z3-Zd))^n;
    B(i+1)=4*km/(1+km)*(gama*(-1)^n/n*(R3/(Z3-Zd))^n);
    i=i+2;
end

%压电纤维2与基体应力连续条件、位移连续条件，负一次多项式；
j=1;
%系数ak
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R3*k/R1*Cji(k+1,2)*(R1/conj(Z3-Z1))^(k+3)*(R3/R1)^2+Z3*k/R1*Cji(k+1,1)*(R1/conj(Z3-Z1))^(k+2)*R3/R1;
   A(i+1,j)=-4/(1+km)*(-R3*k/R1*Cji(k+1,2)*(R1/conj(Z3-Z1))^(k+3)*(R3/R1)^2+Z3*k/R1*Cji(k+1,1)*(R1/conj(Z3-Z1))^(k+2)*R3/R1);
    j=j+1;

end
%系数bk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-R3*k/R2*Cji(k+1,2)*(R2/conj(Z3-Z2))^(k+3)*(R3/R2)^2+Z3*k/R2*Cji(k+1,1)*(R2/conj(Z3-Z2))^(k+2)*R3/R2;
   A(i+1,j)=-4/(1+km)*(-R3*k/R2*Cji(k+1,2)*(R2/conj(Z3-Z2))^(k+3)*(R3/R2)^2+Z3*k/R2*Cji(k+1,1)*(R2/conj(Z3-Z2))^(k+2)*R3/R2);
    j=j+1;

end

%系数ik
for k=1:1:N
    if k==1
        A(i,j)=1;
       A(i+1,j)=4*km/(1+km);
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数ck
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R1/conj(Z3-Z1))^(k+1)*R3/R1;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R1/conj(Z3-Z1))^(k+1)*R3/R1);
    j=j+1;
end
%系数dk
for k=1:1:N
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=-Cji(k,1)*(R2/conj(Z3-Z2))^(k+1)*R3/R2;
   A(i+1,j)=-4/(1+km)*(-Cji(k,1)*(R2/conj(Z3-Z2))^(k+1)*R3/R2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j=j+1;
end
%系数jk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
end

%系数ek
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
  
end
%系数fk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
%系数gk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
%系数hk
for k=1:1:N
    A(i,j)=0;
   A(i+1,j)=0;
    j=j+1;
    A(i,j)=0; 
   A(i+1,j)=0;
    j=j+1;

end
%系数kk
for k=1:1:N
    if k==3
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-3;
       A(i+1,j)=4/(1+kf)*CE*3;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end
%系数lk
for k=1:1:N
    if k==1
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-1;
       A(i+1,j)=4/(1+kf)*CE*1;
        j=j+1;
    else
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
end

B(i)=-conj(gama)*(R3/conj(Z3-Zd))^3+Z3*conj(gama)/R3*(R3/conj(Z3-Zd))^2-(B2-1i*C2)*R3-gama*R3/conj(Z3-Zd)-conj(gama)*Zd/R3*(R3/conj(Z3-Zd))^2;
B(i+1)=-4/(1+km)*(-conj(gama)*(R3/conj(Z3-Zd))^3+Z3*conj(gama)/R3*(R3/conj(Z3-Zd))^2-(B2-1i*C2)*R3-gama*R3/conj(Z3-Zd)-conj(gama)*Zd/R3*(R3/conj(Z3-Zd))^2);
%压电纤维1与基体应力连续条件、位移连续条件，负n次多项式；
i=i+2;

for n=2:1:N
    j=1;
    %系数ak
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-R3*k/R1*(-1)^(n+1)*Cji(k+1,n+1)*(R1/conj(Z3-Z1))^(k+n+2)*(R3/R1)^(n+1)-Z3*k/R1*(-1)^n*Cji(k+1,n)*(R1/conj(Z3-Z1))^(k+n+1)*(R3/R1)^n;
       A(i+1,j)=-4/(1+km)*(-R3*k/R1*(-1)^(n+1)*Cji(k+1,n+1)*(R1/conj(Z3-Z1))^(k+n+2)*(R3/R1)^(n+1)-Z3*k/R1*(-1)^n*Cji(k+1,n)*(R1/conj(Z3-Z1))^(k+n+1)*(R3/R1)^n);
        j=j+1;
    end
    %系数bk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=-R3*k/R2*(-1)^(n+1)*Cji(k+1,n+1)*(R2/conj(Z3-Z2))^(k+n+2)*(R3/R2)^(n+1)-Z3*k/R2*(-1)^n*Cji(k+1,n)*(R2/conj(Z3-Z2))^(k+n+1)*(R3/R2)^n;
       A(i+1,j)=-4/(1+km)*(-R3*k/R2*(-1)^(n+1)*Cji(k+1,n+1)*(R2/conj(Z3-Z2))^(k+n+2)*(R3/R2)^(n+1)-Z3*k/R2*(-1)^n*Cji(k+1,n)*(R2/conj(Z3-Z2))^(k+n+1)*(R3/R2)^n);
        j=j+1;
    end
    
    %系数ik
    for k=1:1:N
        if k==n
            A(i,j)=1;
            A(i+1,j)=4*km/(1+km)*1;
            j=j+1;
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
        else
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
            A(i+1,j)=0;
            j=j+1;
        end
    end
%系数ck
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R1/conj(Z3-Z1))^(k+n)*(R3/R1)^n;       
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R1/conj(Z3-Z1))^(k+n)*(R3/R1)^n);
        j=j+1;
    end
    %系数dk
    
    for k=1:1:N
        A(i,j)=0; 
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=(-1)^n*Cji(k,n)*(R2/conj(Z3-Z2))^(k+n)*(R3/R2)^n;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
       A(i+1,j)=-4/(1+km)*((-1)^n*Cji(k,n)*(R2/conj(Z3-Z2))^(k+n)*(R3/R2)^n);
        j=j+1;
    end    
    %系数jk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end

    %系数ek
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
    end
    %系数fk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end
     %系数gk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end
    %系数hk
    for k=1:1:N
        A(i,j)=0;
       A(i+1,j)=0;
        j=j+1;
        A(i,j)=0;       
       A(i+1,j)=0;
        j=j+1;
    end   
    %系数kk
    for k=1:1:N
        if k==n+2
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-(n+2);
           A(i+1,j)=4/(1+kf)*CE*(n+2);
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    %系数lk
    for k=1:1:N
        if k==n
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=-1;
           A(i+1,j)=4/(1+kf)*CE*1;
            j=j+1;
        else
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
            A(i,j)=0;
           A(i+1,j)=0;
            j=j+1;
        end
    end
    
    B(i)=-conj(gama)*(-1)^(n+1)*(R3/conj(Z3-Zd))^(n+2)-Z3*conj(gama)*(-1)^n/R3*(R3/conj(Z3-Zd))^(n+1)+gama*(-1)^n/n*(R3/conj(Z3-Zd))^n+conj(gama)*(-1)^n*Zd/R3*(R3/conj(Z3-Zd))^(n+1);
    B(i+1)=-4/(1+km)*(-conj(gama)*(-1)^(n+1)*(R3/conj(Z3-Zd))^(n+2)-Z3*conj(gama)*(-1)^n/R3*(R3/conj(Z3-Zd))^(n+1)+gama*(-1)^n/n*(R3/conj(Z3-Zd))^n+conj(gama)*(-1)^n*Zd/R3*(R3/conj(Z3-Zd))^(n+1));
    i=i+2;
end
