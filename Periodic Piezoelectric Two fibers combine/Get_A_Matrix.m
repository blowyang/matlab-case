function A_Matrix=Get_A_Matrix(N,R1,R2,x1,x2,Gm,km,Gf1,kf1,Gf2,kf2)
A_Matrix=zeros(16*N,16*N);
for j=1:N
    [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_One(N,R1,R2,x1,x2,j);
    A_Matrix(j,:)=real(A_Cof_Pos_j);
    A_Matrix(j+N,:)=imag(A_Cof_Pos_j);
    A_Matrix(j+2*N,:)=real(A_Cof_Neg_j);
    A_Matrix(j+3*N,:)=imag(A_Cof_Neg_j);
    [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_Two(N,R1,R2,x1,x2,Gm,km,Gf1,kf1,j);
    A_Matrix(j+4*N,:)=real(A_Cof_Pos_j);
    A_Matrix(j+5*N,:)=imag(A_Cof_Pos_j);
    A_Matrix(j+6*N,:)=real(A_Cof_Neg_j);
    A_Matrix(j+7*N,:)=imag(A_Cof_Neg_j);
    [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_Three(N,R1,R2,x1,x2,j);
    A_Matrix(j+8*N,:)=real(A_Cof_Pos_j);
    A_Matrix(j+9*N,:)=imag(A_Cof_Pos_j);
    A_Matrix(j+10*N,:)=real(A_Cof_Neg_j);
    A_Matrix(j+11*N,:)=imag(A_Cof_Neg_j);
    [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_Four(N,R1,R2,x1,x2,Gm,km,Gf2,kf2,j);
    A_Matrix(j+12*N,:)=real(A_Cof_Pos_j);
    A_Matrix(j+13*N,:)=imag(A_Cof_Pos_j);
    A_Matrix(j+14*N,:)=real(A_Cof_Neg_j);
    A_Matrix(j+15*N,:)=imag(A_Cof_Neg_j);
end

end