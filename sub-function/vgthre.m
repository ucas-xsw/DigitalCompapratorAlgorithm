function vgre=vgthre(x,VDD)
%     global beta_m VTH ku kd
    kd=[0.5	0.33	0.25	4	3	2	1]';
    ku=[2	3	4	0.25	0.33	0.5	1]';
    beta_m=1.4;
    VTH=0.335;
    mfre=sqrt(beta_m*(ku.'*x)/(kd.'*x));
    vgre=(VDD-VTH+VTH*mfre)/mfre;
end
