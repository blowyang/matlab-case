function [A_Matrix,B_Matrix]=Get_AB_Matrix(N,M,R,a,Gm,Gf,km,kf,x)
A_Matrix=zeros(8*N,8*N);
B_Matrix=zeros(8*N,4*M);
for j=1:N
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_One(N,M,R,a,j,x);
    A_Matrix(j,:)=A_Cof_Re;
    B_Matrix(j,:)=B_Cof_Re;
    A_Matrix(j+N,:)=A_Cof_Im;
    B_Matrix(j+N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_One(N,M,R,a,-j,x);
    A_Matrix(j+2*N,:)=A_Cof_Re;
    B_Matrix(j+2*N,:)=B_Cof_Re;
    A_Matrix(j+3*N,:)=A_Cof_Im;
    B_Matrix(j+3*N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Two(N,M,R,a,Gm,Gf,km,kf,j,x);
    A_Matrix(j+4*N,:)=A_Cof_Re;
    B_Matrix(j+4*N,:)=B_Cof_Re;
    A_Matrix(j+5*N,:)=A_Cof_Im;
    B_Matrix(j+5*N,:)=B_Cof_Im;
    [A_Cof_Re,A_Cof_Im,B_Cof_Re,B_Cof_Im]=Get_AB_Matrix_Cof_Two(N,M,R,a,Gm,Gf,km,kf,-j,x);
    A_Matrix(j+6*N,:)=A_Cof_Re;
    B_Matrix(j+6*N,:)=B_Cof_Re;
    A_Matrix(j+7*N,:)=A_Cof_Im;
    B_Matrix(j+7*N,:)=B_Cof_Im;
end

end