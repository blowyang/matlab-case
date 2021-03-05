function B_Matrix=Get_B_Matrix(N,M,km,R,a)
%lambda1=R1/a 
%gama1=G0/G1
B_Matrix=zeros(8*N,4*M);
%第2个方程
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
%第2个方程
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
%第3个方程
for j=1:N
    for k=1:M
        B_Matrix(j+2*N,2*k-1)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+2*N,2*k-1+2*M)=-Get_P_Coef(k,j,0,R,a);
    end   
end
%第4个方程
for j=1:N
    for k=1:M
        B_Matrix(j+3*N,2*k)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+3*N,2*k+2*M)=-Get_P_Coef(k,j,0,R,a);
    end         
end
%第5个方程
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
%第6个方程
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
%第7个方程
for j=1:N
    for k=1:M
        B_Matrix(j+6*N,2*k-1)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+6*N,2*k-1+2*M)=-Get_P_Coef(k,j,0,R,a);
    end   
end
%第8个方程
for j=1:N
    for k=1:M
        B_Matrix(j+7*N,2*k)=-R*Get_P_Diff_Coef(k,j+1,0,R,a);
    end 
    for k=1:M
        B_Matrix(j+7*N,2*k+2*M)=-Get_P_Coef(k,j,0,R,a);
    end         
end

end