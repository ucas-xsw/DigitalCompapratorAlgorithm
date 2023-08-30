function mu=getmu(x)
global mode gi bb
    if mode=="example"
        mug=getmugi1(gi*x',bb);
    else
        mug=zeros(5,1);
        mug(1)=getmugi(x,tr(x),taomax);
        mug(2)=getmugi(x,tl(x),taomax);
        mug(3)=getmugi(x,pw(x),pmax);
        mug(4)=getmugi(x,area(x),pmax);
        mug(5)=getmugi(x,eta(x),etamax);

    end
    mu=sum(mug);
end