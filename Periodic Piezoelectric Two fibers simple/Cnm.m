function res=Cnm(n,m)
tem_1=1;
tem_2=1;
for i=1:1:m
    tem_1=tem_1*(n-(i-1));
    tem_2=tem_2*i;
end
res=tem_1/tem_2;