function P_Coef=Get_P_Coef(k,j,x,R,a)


xN=2*pi;
x0=0;

hn=(xN-x0);
z0=x+R*exp(1i*x0);
zN=x+R*exp(1i*xN);
Tn=hn/2*(Get_P(k,a,z0)*exp(-j*1i*x0)+Get_P(k,a,zN)*exp(-j*1i*xN));
h2n=hn/2;
xc=x0+h2n;
zc=x+R*exp(1i*xc);
T2n=Tn/2+h2n*Get_P(k,a,zc)*exp(-j*1i*xc);
Sn=4/3*T2n-1/3*Tn;
% disp("Sn is ");
% disp(Sn);
hn=h2n;
Tn=T2n;

for np=1:3
    h2n=hn/2;
    xc=x0+h2n;
    T2n=Tn/2;
    while xc<xN
        zc=x+R*exp(1i*xc);
        T2n=T2n+h2n*Get_P(k,a,zc)*exp(-j*1i*xc);
        xc=xc+hn;
    end
    S2n=4/3*T2n-1/3*Tn;   
    Sn=S2n;
    Tn=T2n;
    hn=h2n;     
end
% disp("Sn is ");
% disp(Sn);
check=true;
con=0;
while check
    con=con+1;
    %disp(con);
    h2n=hn/2;
    xc=x0+h2n;
    T2n=Tn/2;
    while xc<xN
        zc=x+R*exp(1i*xc);
        T2n=T2n+h2n*Get_P(k,a,zc)*exp(-j*1i*xc);
        xc=xc+hn;
    end
    S2n=4/3*T2n-1/3*Tn;
%     disp("S2n is ");
%     disp(S2n);
    if(abs(real(S2n-Sn)/real(Sn))<0.001 || abs(real(S2n))<1e-10)
        check=false;
        P_Coef=real(S2n);
    else
        Sn=S2n;
        Tn=T2n;
        hn=h2n;
    end       
end
P_Coef=P_Coef/2/pi;

% NP=500;
% h=2*pi/NP;
% zd=x;
% inte=0;
% z0=zd+R*exp(1i*0);
% inte=inte+h/6*Get_P(k,a,z0)*exp(-j*1i*0);
% zNP=zd+R*exp(1i*2*pi);
% inte=inte+h/6*Get_P(k,a,zNP)*exp(-j*1i*2*pi);
% 
% for np=1:(NP-1)
%     zi=zd+R*exp(1i*np*h);
%     inte=inte+h/3*Get_P(k,a,zi)*exp(-j*1i*np*h);
% end
% 
% for np=1:NP
%     zi=zd+R*exp(1i*(np-0.5)*h);
%     inte=inte+2*h/3*Get_P(k,a,zi)*exp(-j*1i*(np-0.5)*h);
% end
% P_Coef=real(inte)/2/pi;

end