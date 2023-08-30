function mugi=getmugi1(gx,b)
global b3
   if  size(gx,2)>1
       mugi=zeros(size(gx,1),size(gx,2));
       b2=b*ones(1,size(gx,2));
       h=(gx>=0) & (gx<=b2);
       mugi(h)=1-gx(h)./b2(h);  
       for i=1:size(gx,2)
           for j=1:size(gx,1)
               gij=gx(j,i)-b2(j,i);
               mugi(j,i);
               if gij<=0
                   mugi(j,i)=1;
               else
                   mugi(j,i)=0;
               end
%                if gij>=0 && gij<=b3
%                    mugi(j,i)=1-gij/b3;
%                end
           end
       end
   else
       mugi=zeros(size(gx,1),1);
       mugi(gx<=b)=1;
%        h=(gx>=0) & (gx<=b);
%        mugi(h)=1-gx(h)./b(h);
   end
end