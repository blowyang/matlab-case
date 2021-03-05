function [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,j)

NP=100;
h=2*pi/NP;
zd=x2;
inte_A_Cof=0;
inte_B_Cof=0;
z0=zd+R2*exp(1i*0);
[A_Cof,B_Cof]=AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,z0);
inte_A_Cof=inte_A_Cof+h/6*A_Cof*exp(-j*1i*0);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*0);
zNP=zd+R2*exp(1i*2*pi);
[A_Cof,B_Cof]=AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,zNP);

inte_A_Cof=inte_A_Cof+h/6*A_Cof*exp(-j*1i*2*pi);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*2*pi);

for np=1:(NP-1)
    zi=zd+R2*exp(1i*np*h);
    [A_Cof,B_Cof]=AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,zi);
    inte_A_Cof=inte_A_Cof+h/3*A_Cof*exp(-j*1i*np*h);
    inte_B_Cof=inte_B_Cof+h/3*B_Cof*exp(-j*1i*np*h);
end

for np=1:NP
    zi=zd+R2*exp(1i*(np-0.5)*h);
    [A_Cof,B_Cof]=AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,zi);
    inte_A_Cof=inte_A_Cof+2*h/3*A_Cof*exp(-j*1i*(np-0.5)*h);
    inte_B_Cof=inte_B_Cof+2*h/3*B_Cof*exp(-j*1i*(np-0.5)*h);
end
A_Cof_Re=real(inte_A_Cof)/2/pi;
B_Cof_Re=real(inte_B_Cof)/2/pi;
A_Cof_Im=imag(inte_A_Cof)/2/pi;
B_Cof_Im=imag(inte_B_Cof)/2/pi;


% [r,c]=size(A_Cof_Re);
% for i=1:r
%     for j=1:c
%         temp=A_Cof_Re(i,j);
%         if abs(temp)<1e-10
%             A_Cof_Re(i,j)=0;
%         end
%     end
% end
% [r,c]=size(B_Cof_Re);
% for i=1:r
%     for j=1:c
%         temp=B_Cof_Re(i,j);
%         if abs(temp)<1e-10
%             B_Cof_Re(i,j)=0;
%         end
%     end
% end
% [r,c]=size(A_Cof_Im);
% for i=1:r
%     for j=1:c
%         temp=A_Cof_Im(i,j);
%         if abs(temp)<1e-10
%             A_Cof_Im(i,j)=0;
%         end
%     end
% end
% [r,c]=size(B_Cof_Im);
% for i=1:r
%     for j=1:c
%         temp=B_Cof_Im(i,j);
%         if abs(temp)<1e-10
%             B_Cof_Im(i,j)=0;
%         end
%     end
% end
%%
% xN=2*pi;
% x0=0;
% 
% hn=(xN-x0);
% z0=x+R*exp(1i*x0);
% zN=x+R*exp(1i*xN);
% Tn=hn/2*(A_Matrix_Cof(R,z0)*exp(-j*1i*x0)+A_Matrix_Cof(R,zN)*exp(-j*1i*xN));
% h2n=hn/2;
% xc=x0+h2n;
% zc=x+R*exp(1i*xc);
% T2n=Tn/2+h2n*A_Matrix_Cof(R,zc)*exp(-j*1i*xc);
% Sn=4/3*T2n-1/3*Tn;
% % disp("Sn is ");
% % disp(Sn);
% hn=h2n;
% Tn=T2n;
% 
% for np=1:3
%     h2n=hn/2;
%     xc=x0+h2n;
%     T2n=Tn/2;
%     while xc<xN
%         zc=x+R*exp(1i*xc);
%         T2n=T2n+h2n*A_Matrix_Cof(R,zc)*exp(-j*1i*xc);
%         xc=xc+hn;
%     end
%     S2n=4/3*T2n-1/3*Tn;   
%     Sn=S2n;
%     Tn=T2n;
%     hn=h2n;     
% end
% % disp("Sn is ");
% % disp(Sn);
% check=true;
% %con=0;
% while check
%     %con=con+1;
%     %disp(con);
%     h2n=hn/2;
%     xc=x0+h2n;
%     T2n=Tn/2;
%     while xc<xN
%         zc=x+R*exp(1i*xc);
%         T2n=T2n+h2n*A_Matrix_Cof(R,zc)*exp(-j*1i*xc);
%         xc=xc+hn;
%     end
%     S2n=4/3*T2n-1/3*Tn;
% %     disp("S2n is ");
% %     disp(S2n);
%     if(abs(real(S2n-Sn)/real(Sn))<0.001 )
%         check=false;
%         A_Cof=real(S2n);
%     elseif(abs(real(S2n))<1e-10)
%         check=false;
%         A_Cof=0;
%     else
%         Sn=S2n;
%         Tn=T2n;
%         hn=h2n;
%     end       
% end
% A_Cof=A_Cof/2/pi;
%%
% NP=100;
% h=2*pi/NP;
% zd=x;
% inte=0;
% z0=zd+R*exp(1i*0);
% inte=inte+h/6*A_Matrix_Cof(N,M,R,a,z0)*exp(-j*1i*0);
% zNP=zd+R*exp(1i*2*pi);
% inte=inte+h/6*A_Matrix_Cof(N,M,R,a,zNP)*exp(-j*1i*2*pi);
% 
% for np=1:(NP-1)
%     zi=zd+R*exp(1i*np*h);
%     inte=inte+h/3*A_Matrix_Cof(N,M,R,a,zi)*exp(-j*1i*np*h);
% end
% 
% for np=1:NP
%     zi=zd+R*exp(1i*(np-0.5)*h);
%     inte=inte+2*h/3*A_Matrix_Cof(N,M,R,a,zi)*exp(-j*1i*(np-0.5)*h);
% end
% A_Cof=real(inte)/2/pi;
end