function res=Cji(x,y)
tem_1=1;
tem_2=1;
for i=1:1:y
    tem_1=tem_1*(x+i-1);
    tem_2=tem_2*i;
end
res=tem_1/tem_2;