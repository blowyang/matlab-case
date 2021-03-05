function B_Cof=Get_B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,j)

NP=80;
h=2*pi/NP;
zd=x1;
inte_B_Cof=0;
z0=zd+R1*exp(1i*0);
B_Cof=B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,z0);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*0);
zNP=zd+R1*exp(1i*2*pi);
B_Cof=B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,zNP);
inte_B_Cof=inte_B_Cof+h/6*B_Cof*exp(-j*1i*2*pi);

for np=1:(NP-1)
    zi=zd+R1*exp(1i*np*h);
    B_Cof=B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,zi);
    inte_B_Cof=inte_B_Cof+h/3*B_Cof*exp(-j*1i*np*h);
end

for np=1:NP
    zi=zd+R1*exp(1i*(np-0.5)*h);
    B_Cof=B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,zi);
    inte_B_Cof=inte_B_Cof+2*h/3*B_Cof*exp(-j*1i*(np-0.5)*h);
end
B_Cof=inte_B_Cof/2/pi;

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