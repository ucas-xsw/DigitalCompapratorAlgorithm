%% b value test: etamax, taomax, pmax,amax
clear
clc

global VDD VTH alp_p  fclk ceq gnd ag  gpu beta_m ku kd vgthid
VDD=1.2;
alp_p=0.5;
fclk=500e6;
beta_m=1.4;
m=3;
H=[1.40E+00	1.40E+00	1.40E+00	1.40E+00	1.40E+00	1.40E+00	1.40E+00]'*1e-6;
W=[8.00E-01	1.00E+00	1.40E+00	8.00E-01	1.00E+00	1.40E+00	6.00E-01]'*1e-6;;
trs=[2.86E-02	4.96E-02	7.45E-02	1.71E-02	1.96E-02	2.20E-02	1.39E-02]'*1e-9;
tfs=[1.44E-02	1.66E-02	1.80E-02	2.33E-02	3.75E-02	5.47E-02	1.24E-02]'*1e-9;
Ih=[1.79E+02	2.69E+02	3.59E+02	5.09E+01	3.57E+01	2.80E+01	8.98E+01]'*1e-6;
Il=[2.08E+01	3.12E+01	4.16E+01	2.71E+00	1.50E+00	1.06E+00	1.04E+01]'*1e-6;
tp=(trs+tfs)/2;
req=VDD./(Il+Ih);
geq=(Il+Ih)/VDD;
ag=H.*W;
cin=tp./req;
geqp=cin./trs;
geqn=cin./tfs;
ceq=cin;
gnd=geqn;
gpu=geqp;
kd=[0.5	0.33	0.25	4	3	2	1]';
ku=[2	3	4	0.25	0.33	0.5	1]';

VTH=0.335;
etamax=0.1;
taomax=0.2e-9;
vgthid=0.5;
eta1=((VDD-VTH)/(sqrt(etamax)+vgthid-VTH))^2*kd'-beta_m*ku';
eta2=-((VDD-VTH)/(-sqrt(etamax)+vgthid-VTH))^2*kd'+beta_m*ku';
tao1=0.345*ceq'-taomax*geqp';
tao2=0.345*ceq'-taomax*geqn';

x=1*ones(2*m+1,1);
pmin=pw(x)*1e6 % uw
trxmin=tr(x)*1e9%ns
tlxmin=tl(x)*1e9%ns
etaxmin=eta(x)

axmin=area(x)*1e12%um2

x=10*ones(2*m+1,1);
pmax=pw(x)*1e6 % uw
trxmax=tr(x)*1e9%ns
tlxmax=tl(x)*1e9%ns
etaxmax=eta(x)
axmax=area(x)*1e12%um2
Mg=7;
k=0;
Mf=zeros(1e7,1);
vgthid_vec=zeros(1e7,1);
% for i1=1:Mg
%     for i2=1:Mg
%         for i3=1:Mg
%             for i4=1:Mg
%                     for i5=1:Mg
%                         for i6=1:Mg
%                              for i7=1:Mg
%                                 x(1)=i1;
%                                 x(2)=i2;
%                                 x(3)=i3;
%                                 x(4)=i4;     
%                                 x(5)=i5;
%                                 x(6)=i6;
%                                 x(7)=i7;
%                                 k=k+1;
%                                 Mf(k)=sqrt(beta_m*ku'*x/(kd'*x));
%                                 vgthid_vec(k)= (VDD-VTH+Mf(k)*VTH)/Mf(k);
%                              end
%                         end    
%                     end        
%             end
%         end
%     end
% end
% vgthre=sort(vgthid_vec,'descend');
% vgthrere_1=vgthre((vgthre<1.2)&(vgthre>0));
beta_m=10%beta_m<24
Mg=4;
x=[Mg Mg Mg 0 0 0 0]';
Mf1=sqrt(beta_m*ku'*x/(kd'*x))
vgthid_1= (VDD-VTH+Mf1*VTH)/Mf1
x=[0 0 0 Mg Mg Mg 0]';
Mf2=sqrt(beta_m*ku'*x/(kd'*x))
vgthid_2= (VDD-VTH+Mf2*VTH)/Mf2
vgthid=0.5;
Mf3=(VDD-VTH)/(vgthid-VTH)


















