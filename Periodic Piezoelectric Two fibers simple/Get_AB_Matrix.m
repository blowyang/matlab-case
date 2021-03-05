function [A_Matrix,B_Matrix]=Get_AB_Matrix(N,M,R1,R2,x1,x2,a,Gm,km,Gf1,kf1,Gf2,kf2)
A_Matrix=zeros(16*N,16*N);
B_Matrix=zeros(16*N,4*M);
for j=1:N
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,j);
    A_Matrix(j,:)=A_Cof_Re;
    B_Matrix(j,:)=B_Cof_Re;
    A_Matrix(j+N,:)=A_Cof_Im;
    B_Matrix(j+N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,-j);
    A_Matrix(j+2*N,:)=A_Cof_Re;
    B_Matrix(j+2*N,:)=B_Cof_Re;
    A_Matrix(j+3*N,:)=A_Cof_Im;
    B_Matrix(j+3*N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,Gm,Gf1,km,kf1,j);
    A_Matrix(j+4*N,:)=A_Cof_Re;
    B_Matrix(j+4*N,:)=B_Cof_Re;
    A_Matrix(j+5*N,:)=A_Cof_Im;
    B_Matrix(j+5*N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,Gm,Gf1,km,kf1,-j);
    A_Matrix(j+6*N,:)=A_Cof_Re;
    B_Matrix(j+6*N,:)=B_Cof_Re;
    A_Matrix(j+7*N,:)=A_Cof_Im;
    B_Matrix(j+7*N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,j);
    A_Matrix(j+8*N,:)=A_Cof_Re;
    B_Matrix(j+8*N,:)=B_Cof_Re;
    A_Matrix(j+9*N,:)=A_Cof_Im;
    B_Matrix(j+9*N,:)=B_Cof_Im;    
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,-j);
    A_Matrix(j+10*N,:)=A_Cof_Re;
    B_Matrix(j+10*N,:)=B_Cof_Re;
    A_Matrix(j+11*N,:)=A_Cof_Im;
    B_Matrix(j+11*N,:)=B_Cof_Im;     
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,j);
    A_Matrix(j+12*N,:)=A_Cof_Re;
    B_Matrix(j+12*N,:)=B_Cof_Re;
    A_Matrix(j+13*N,:)=A_Cof_Im;
    B_Matrix(j+13*N,:)=B_Cof_Im;      
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,Gm,Gf2,km,kf2,-j);
    A_Matrix(j+14*N,:)=A_Cof_Re;
    B_Matrix(j+14*N,:)=B_Cof_Re;
    A_Matrix(j+15*N,:)=A_Cof_Im;
    B_Matrix(j+15*N,:)=B_Cof_Im;            
end

end