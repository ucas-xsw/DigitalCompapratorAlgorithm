function select=alphaleq(x,u) % x的第1维不是确定的
global alpha
    select=2;
    mu1=getmu(x);
    mu2=getmu(u);
    f1=EDP(x);
    f2=EDP(u);    
    if mu1>mu2
        select=1;
     end
     if (mu1>=alpha) & (mu2>=alpha)
            if f1<f2
                select=1;
            else
                select=2;
            end   
     end
     if mu1==mu2
            if f1<f2
                select=1;
            else
                select=2;
            end 
     end
end