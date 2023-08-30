function  et=geteta(x)  
    global vgthid VDD
        et=(vgthre(x,VDD)-vgthid)^2;
    
end