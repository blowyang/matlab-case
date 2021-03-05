function A_Matrix=Get_A_Matrix(N,Gm,Gf,kf,km)
A_Matrix=zeros(8*N,8*N);

for j=1:N
    if(j<3)
        if(j==1)
            A_Matrix(j,1)=2;
            A_Matrix(j,1+6*N)=-1;
        end
        if(j==2)
            A_Matrix(j,3)=1;            
            A_Matrix(j,3+6*N)=-1;
        end        
    else
        A_Matrix(j,2*j-1)=1;
        A_Matrix(j,2*j-5+4*N)=(j-2);
        A_Matrix(j,2*j-1+6*N)=-1;
    end
end

for j=1:N
    if(j<3)
        if(j==1)
            A_Matrix(j+N,2+6*N)=-1;
        end
        if(j==2)
            A_Matrix(j+N,4)=1;
            A_Matrix(j+N,4+6*N)=1;
        end
    else
        A_Matrix(j+N,2*j)=1;
        A_Matrix(j+N,2*j-4+4*N)=-(j-2);
        A_Matrix(j+N,2*j+6*N)=1;
    end  
end

for j=1:N
    if(2*j+3<=2*N)
        A_Matrix(j+2*N,2*j+3)=(j+2);
    end
    A_Matrix(j+2*N,2*j-1+2*N)=1;
    A_Matrix(j+2*N,2*j-1+4*N)=-1;
end
%第4个方程
for j=1:N
    if(2*j+4<=2*N)
        A_Matrix(j+3*N,2*j+4)=(j+2);
    end
    A_Matrix(j+3*N,2*j+2*N)=1;
    A_Matrix(j+3*N,2*j+4*N)=1;
end
%第5个方程
for j=1:N
    if(j<3)
        if(j==1)
            A_Matrix(j+4*N,1)=Gm/Gf*(kf-1);
            A_Matrix(j+4*N,1+6*N)=1;
        end
        if(j==2)
            A_Matrix(j+4*N,3)=Gm/Gf*kf;
            A_Matrix(j+4*N,3+6*N)=1;
        end        
    else
        A_Matrix(j+4*N,2*j-1)=Gm/Gf*kf;
        A_Matrix(j+4*N,2*j-5+4*N)=-(j-2);
        A_Matrix(j+4*N,2*j-1+6*N)=1;
    end
end
%第6个方程
for j=1:N
    if(j<3)
        if(j==1)
            A_Matrix(j+5*N,2)=Gm/Gf*(kf+1);
            A_Matrix(j+5*N,2+6*N)=-1;
        end
        if(j==2)
            A_Matrix(j+5*N,4)=Gm/Gf*kf;
            A_Matrix(j+5*N,4+6*N)=-1;
        end
    else
        A_Matrix(j+5*N,2*j)=Gm/Gf*kf;
        A_Matrix(j+5*N,2*j-4+4*N)=(j-2);
        A_Matrix(j+5*N,2*j+6*N)=-1;
    end  
end
%第7个方程
for j=1:N
    if(2*j+3<=2*N)
        A_Matrix(j+6*N,2*j+3)=Gm/Gf*(j+2);
    end
    A_Matrix(j+6*N,2*j-1+2*N)=Gm/Gf;
    A_Matrix(j+6*N,2*j-1+4*N)=km;
end
%第8个方程
for j=1:N
    if(2*j+4<=2*N)
        A_Matrix(j+7*N,2*j+4)=Gm/Gf*(j+2);
    end
    A_Matrix(j+7*N,2*j+2*N)=Gm/Gf;
    A_Matrix(j+7*N,2*j+4*N)=-km;
end
end