function B_Matrix=Get_B_Matrix(N,M,R1,R2,x1,x2,a,km)
B_Matrix=zeros(16*N,4*M);
for j=1:N
    B_Cof=Get_B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,j);
    B_Matrix(j,:)=real(B_Cof);
    B_Matrix(j+N,:)=imag(B_Cof);  
    B_Cof=Get_B_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,-j);
    B_Matrix(j+2*N,:)=real(B_Cof);
    B_Matrix(j+3*N,:)=imag(B_Cof);
    
    B_Cof=Get_B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,j);
    B_Matrix(j+4*N,:)=real(B_Cof);
    B_Matrix(j+5*N,:)=imag(B_Cof);
    B_Cof=Get_B_Matrix_Cof_Two(N,M,R1,R2,x1,x2,a,km,-j);
    B_Matrix(j+6*N,:)=real(B_Cof);
    B_Matrix(j+7*N,:)=imag(B_Cof);
    
    B_Cof=Get_B_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,j);
    B_Matrix(j+8*N,:)=real(B_Cof);
    B_Matrix(j+9*N,:)=imag(B_Cof);
    B_Cof=Get_B_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,-j);
    B_Matrix(j+10*N,:)=real(B_Cof);
    B_Matrix(j+11*N,:)=imag(B_Cof);
    
    B_Cof=Get_B_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,km,j);
    B_Matrix(j+12*N,:)=real(B_Cof);
    B_Matrix(j+13*N,:)=imag(B_Cof); 
    B_Cof=Get_B_Matrix_Cof_Four(N,M,R1,R2,x1,x2,a,km,-j);
    B_Matrix(j+14*N,:)=real(B_Cof);
    B_Matrix(j+15*N,:)=imag(B_Cof);         
end

end