function P=Get_P(n,a,z)
m3=-1/6;
m7=1/56;
coef=143/(84*a);
%coef=1;
ca=0;
if(n==1)
    ca=1;
end
if (n>1)&&(n<4)
    ca=2;
end
if(n==4)
    ca=3;
end
if(n>4)&&(n<8)
    ca=4;
end
if(n==8)
    ca=5;
end
if(n>8)
    ca=6;
end
switch ca
    case 1
        P=coef*z;
    case 2
        P=Get_P(1,a,z)*Get_P(n-1,a,z);
    case 3
        P=Get_P(1,a,z)*Get_P(n-1,a,z)-n*m3;
    case 4
        P=Get_P(1,a,z)*Get_P(n-1,a,z)-m3*Get_P(n-4,a,z);
    case 5
        P=Get_P(1,a,z)*Get_P(n-1,a,z)-m3*Get_P(n-4,a,z)-n*m7;
    case 6
        P=Get_P(1,a,z)*Get_P(n-1,a,z)-m3*Get_P(n-4,a,z)-m7*Get_P(n-8,a,z);
end
end