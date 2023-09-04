%% point: D*M

function xbest=discrete_nelder_mead(point)
global ka kb kg D eps M
delta=1;
xbest=zeros(1,M);
v=point;
u=zeros(D+1,M);
i=0;
while delta>eps
    edpp=EDP(point);
    edph(1,:)=max(edpp);
    edpl=min(edpp);
    l=(edpp==edpl);
    h=find(edpp==edph(1,:));
    ph=point(h,:);
    pl=point(l,:);
    pavg=sum(point)/D-sum(point(h))/D;
    ku=floor(norm(pavg-ph(1,:),2));
    pref=ph(1,:)+ka*ku*sign(pavg-ph(1,:));
    select=alphaleq(pref,pl);
    u=xbest;
   if select==1
       pexp=pref+kb*ku*sign(pavg-ph(1,:));
        if alphaleq(pexp,pl)==1
            ph=ones(size(ph,1),1)*pexp;
        else
            ph=ones(size(ph,1),1)*pref;
        end
   else
%         point(1:h-1,:);
%         point(h+1:end,:);
       pnh=point;
       pnh(h,:)=[];
       if alphaleq(pnh,pref)==1
           if alphaleq(ph(1,:),pref)==1
                pcon=pref-kg*ku*sign(pavg-ph(1,:));
           else
                ph=ones(size(ph,1),1)*pref;
                pcon=pref-kg*ku*sign(pavg-ph(1,:));
           end
           
           if alphaleq(ph(1,:),pcon)==1
               point=round((point+ones(D+1,1)*pl(1,:))/2);
           else
               ph=ones(size(ph,1),1)*pcon;
           end
       else
           ph=ones(size(ph,1),1)*pref;
       end
   end
    i=i+1
   l=(EDP(point))==min(EDP(point));
   pl=point(l,:);
   point(h,:)=ph;
   point(l,:)=pl;
%    v=[pavg,point']';
   v=point;
   vl=v((EDP(v))==min(EDP(v)),:);
   xbest=vl(1,:);
   min(EDP(v))
   delta=norm(EDP(u)-EDP(xbest),2)
   %delta=sqrt((EDP(v)-mean (EDP(v)))'*(EDP(v)-mean (EDP(v)))/(D+1))
end
