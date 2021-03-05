function B_Matrix=Get_B_Matrix(N,M,km,R,a)
%lambda1=R1/a 
%gama1=G0/G1
B_Matrix=zeros(8*N,4*M);
%��2������
for j=1:N
    if(j==1)
        for k=1:M
            B_Matrix(j,2*k-1)=-(Get_P_Coef(k,1,0,R,a)+R*Get_P_Diff_Coef(k,0,0,R,a));
        end       
    else
        for k=1:M
            B_Matrix(j,2*k-1)=-Get_P_Coef(k,j,0,R,a);
        end  
    end  
end
%��2������
for j=1:N
    if(j==1)
        for k=1:M
            B_Matrix(j+N,2*k)=(Get_P_Coef(k,1,0,R,a)-R*Get_P_Diff_Coef(k,0,0,R,a));
        end
    else
        for k=1:M
            B_Matrix(j+N,2*k)=-Get_P_Coef(k,j,0,R,a);
        end
    end  
end
%��3������
for j=1:N
    for k=1:M
        B_Matrix(j+2*N,2*k-1)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+2*N,2*k-1+2*M)=-Get_P_Coef(k,j,0,R,a);
    end   
end
%��4������
for j=1:N
    for k=1:M
        B_Matrix(j+3*N,2*k)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+3*N,2*k+2*M)=-Get_P_Coef(k,j,0,R,a);
    end         
end
%��5������
for j=1:N
    if(j==1)
        for k=1:M
            B_Matrix(j+4*N,2*k-1)=-(km*Get_P_Coef(k,1,0,R,a)-R*Get_P_Diff_Coef(k,0,0,R,a));
        end       
    else
        for k=1:M
            B_Matrix(j+4*N,2*k-1)=-km*Get_P_Coef(k,j,0,R,a);
        end  
    end  
end
%��6������
for j=1:N
    if(j==1)
        for k=1:M
            B_Matrix(j+5*N,2*k)=-(km*Get_P_Coef(k,1,0,R,a)+R*Get_P_Diff_Coef(k,0,0,R,a));
        end       
    else
        for k=1:M
            B_Matrix(j+5*N,2*k)=-km*Get_P_Coef(k,j,0,R,a);
        end  
    end  
end
%��7������
for j=1:N
    for k=1:M
        B_Matrix(j+6*N,2*k-1)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+6*N,2*k-1+2*M)=-Get_P_Coef(k,j,0,R,a);
    end   
end
%��8������
for j=1:N
    for k=1:M
        B_Matrix(j+7*N,2*k)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+7*N,2*k+2*M)=-Get_P_Coef(k,j,0,R,a);
    end         
end

end