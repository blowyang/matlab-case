function [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
    Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012)

delta_cof=zeros(4,4);

NP=100;
h=a/NP;
%%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,z0);
inte=inte+h/6*Sm22;
x2N=a/2;
zN=x2N+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zN);
inte=inte+h/6*Sm22;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+h/3*Sm22;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+2*h/3*Sm22;
end
delta_cof(1,1)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,z0);
inte=inte+h/6*Sm22;
x2N=a/2;
zN=x2N+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zN);
inte=inte+h/6*Sm22;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+h/3*Sm22;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+2*h/3*Sm22;
end
delta_cof(1,2)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,z0);
inte=inte+h/6*Sm22;
x2N=a/2;
zN=x2N+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zN);
inte=inte+h/6*Sm22;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+h/3*Sm22;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+2*h/3*Sm22;
end
delta_cof(1,3)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,z0);
inte=inte+h/6*Sm22;
x2N=a/2;
zN=x2N+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zN);
inte=inte+h/6*Sm22;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+h/3*Sm22;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+2*h/3*Sm22;
end
delta_cof(1,4)=inte;

%%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(2,1)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(2,2)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(2,3)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(2,4)=inte;
%%
inte=0;
x20=-a/2;
z0=a/2+x20*1i;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,z0);
inte=inte+h/6*Sm11;
x2N=a/2;
zN=a/2+x2N*1i;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zN);
inte=inte+h/6*Sm11;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+h/3*Sm11;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+2*h/3*Sm11;
end
delta_cof(3,1)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,z0);
inte=inte+h/6*Sm11;
x2N=a/2;
zN=x2N*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zN);
inte=inte+h/6*Sm11;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+h/3*Sm11;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+2*h/3*Sm11;
end
delta_cof(3,2)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,z0);
inte=inte+h/6*Sm11;
x2N=a/2;
zN=x2N*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zN);
inte=inte+h/6*Sm11;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+h/3*Sm11;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+2*h/3*Sm11;
end
delta_cof(3,3)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,z0);
inte=inte+h/6*Sm11;
x2N=a/2;
zN=x2N*1i+a/2;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zN);
inte=inte+h/6*Sm11;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+h/3*Sm11;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+2*h/3*Sm11;
end
delta_cof(3,4)=inte;

%%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce1,Xdf1,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(4,1)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce2,Xdf2,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(4,2)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce3,Xdf3,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(4,3)=inte;
%
%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce4,Xdf4,zi);
    inte=inte+2*h/3*Sm12;
end
delta_cof(4,4)=inte;
%%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,z0);
inte=inte+h/6*Sm22;
x2N=a/2;
zN=x2N+a/2*1i;
[Sm22,~,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zN);
inte=inte+h/6*Sm22;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+h/3*Sm22;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+2*h/3*Sm22;
end
Stress_res(1,1)=S022*a-inte;
%%
inte=0;
x20=-a/2;
z0=x20+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N+a/2*1i;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n+a/2*1i;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+2*h/3*Sm12;
end
Stress_res(2,1)=S012*a-inte;
%%
inte=0;
x20=-a/2;
z0=a/2+x20*1i;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,z0);
inte=inte+h/6*Sm11;
x2N=a/2;
zN=a/2+x2N*1i;
[~,Sm11,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zN);
inte=inte+h/6*Sm11;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+h/3*Sm11;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,Sm11,~]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+2*h/3*Sm11;
end
Stress_res(3,1)=S011*a-inte;
%%
inte=0;
x20=-a/2;
z0=x20*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,z0);
inte=inte+h/6*Sm12;
x2N=a/2;
zN=x2N*1i+a/2;
[~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zN);
inte=inte+h/6*Sm12;

for np=1:(NP-1)
    x2n=-a/2+np*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+h/3*Sm12;
end

for np=1:NP
    x2n=-a/2+(np-0.5)*h;
    zi=x2n*1i+a/2;
    [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce0,Xdf0,zi);
    inte=inte+2*h/3*Sm12;
end
Stress_res(4,1)=S012*a-inte;
%%
res=delta_cof\Stress_res;
delta1=res(1,1);
delta2=res(2,1);
delta3=res(3,1);
delta4=res(4,1);
end
