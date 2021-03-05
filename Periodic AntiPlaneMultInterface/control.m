function y=control(x1,x2)
m3=-1/6;
m7=1/56;
% if x==1
%     y=1;
% else
%     y=x*control(x-1);
% end
% switch x
%     case 1
%         y=1;
%     otherwise
%         y=x*control(x-1);
% end
% case fen lei bianliang n
n=2;
if (x2>x1) 
    n=0;   
end
if (x2==x1)
    n=1;
end
% 3 cases
switch n
    case 0
        y=0;
    case 1
        y=1;
    otherwise
        m=0;
        if (x1==1)
            m=1;
        end
        if (x1>1)&&(x1<4)
            m=2;
        end
        if(x1==4)
            m=3;
        end
        if(x1>4)&&(x1<8)
           m=4; 
        end
        if(x1==8)
            m=5;
        end
        switch m
            case 1
                y=0;
            case 2
                if(x2>0)
                    y=control(x1-1,x2-1);
                else
                    y=0;
                end
            case 3
                if(x2>0)
                    y=control(x1-1,x2-1);
                else
                    y=-x1*m3;
                end 
            case 4
                if(x2>0)
                    y=control(x1-1,x2-1)-m3*control(x1-4,x2);
                else
                    y=-m3*control(x1-4,x2);
                end 
            case 5
                if(x2>0)
                    y=control(x1-1,x2-1)-m3*control(x1-4,x2);
                else
                    y=-m3*control(x1-4,x2)-x1*m7;
                end
            otherwise
                if(x2>0)
                    y=control(x1-1,x2-1)-m3*control(x1-4,x2)-m7*control(x1-8,x2);
                else
                    y=-m3*control(x1-4,x2)-m7*control(x1-8,x2);
                end                
        end
end
end
