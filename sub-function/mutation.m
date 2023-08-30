function L=mutation(x)
global NP M alp beta kk
L=zeros(NP,M);
edps=sort(EDP(x));
for i=1:NP
    b(i,:)=(x(EDP(x)==edps(1),:));
end

for i=1:NP
    if(i<=kk)
        a=zeros(i+kk-1,M);
        a(1:i-1,:)=b(1:i-1,:);
        a(i:i+kk-1,:)=b(i+1:i+kk,:);
        
    end
    if(i>kk && i<=NP-kk)
        a=zeros(2*kk,M);
        a(1:kk,:)=b(i-kk:i-1,:);
        a(kk+1:2*kk,:)=b(i+1:i+kk,:);
    end
    if (i>NP-kk)
        a=zeros(kk+NP-i,M);
        a(1:kk,:)=b(i-kk:i-1,:);
        a(kk+1:kk+NP-i,:)=b(i+1:NP,:);     
    end
    xbest=a(EDP(a)==min(EDP(a)),:);
    if(size(a,1)>1)
        num=1:size(a,1);
        xr= a(randperm(numel(num),2),:);
        L(i,:)=x(i,:)+alp*(xbest(1)-x(i,:))+beta*(xr(1,:)-xr(2,:));
    else
        L(i,:)=x(i,:)+alp*(xbest(1)-x(i,:));
    end
end
end