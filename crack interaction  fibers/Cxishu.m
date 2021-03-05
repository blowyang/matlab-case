function [CA, CB]=Cxishu(A ,B,N)
for k=1:1:12*N
    for j=1:1:12*N
        tem_1=A(k,2*j-1);
        tem_2=A(k,2*j);
        CA(k,2*j-1)=tem_1;
        CA(k,2*j)=tem_2;
        CA(12*N+k,2*j-1)=conj(tem_2);
        CA(12*N+k,2*j)=conj(tem_1);
        tem_1=B(k);
        CB(k)=tem_1;
        CB(12*N+k)=conj(tem_1);
    end
end
      