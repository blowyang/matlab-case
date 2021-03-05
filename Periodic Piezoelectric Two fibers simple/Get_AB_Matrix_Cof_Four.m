function [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,j)

NP=100;
h=2*pi/NP;
zd=x2;
inte_A_Cof=0;
inte_B_Cof=0;
z0=zd+R2*exp(1i*0);
[A_Cof,B_Cof]=AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,z0);
inte_A_Cof=inte_A_Cof+h/6*A_Cof*exp(-j*1i*0);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*0);
zNP=zd+R2*exp(1i*2*pi);
[A_Cof,B_Cof]=AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,zNP);
inte_A_Cof=inte_A_Cof+h/6*A_Cof*exp(-j*1i*2*pi);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*2*pi);

for np=1:(NP-1)
    zi=zd+R2*exp(1i*np*h);
    [A_Cof,B_Cof]=AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,zi);
    inte_A_Cof=inte_A_Cof+h/3*A_Cof*exp(-j*1i*np*h);
    inte_B_Cof=inte_B_Cof+h/3*B_Cof*exp(-j*1i*np*h);
end

for np=1:NP
    zi=zd+R2*exp(1i*(np-0.5)*h);
    [A_Cof,B_Cof]=AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,zi);
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

end