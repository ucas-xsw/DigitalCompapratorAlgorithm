function mugi=getmugi(x,gx,b)
   mugi=0;
   if gx<=0
       mugi=1;
   end
   if (gx>=0) & (gx<=b)
       mugi=1-gx/b;       
   end
end