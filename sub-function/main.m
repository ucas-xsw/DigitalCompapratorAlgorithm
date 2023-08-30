
global alp beta xmax xmin NP cr D M N eps ka kb kg mode gi alpha bb b3 kk
%% ADC








%% DENM
M=5; %variable
NP=100;
xmax=99;
xmin=0;
D=100; % DNM
alp =0.8; % mutation
beta =0.5;
cr=0.3; % crossover
ka =1;%  integer reflection coefficien
kb =6;%  Ó°Ïìfmin´óĞ¡£¬integer expansion coefficient
kg =1;%   integer contraction coefficient
xbest_pst=zeros(M,1);
xbest=zeros(M,1); 

mu=0;
sigma=0.5;
N=10;
eps=0.1;
g=2;
xbest(g-1,1:M)=ones(1,M)*100;
mode="example"; % 'objective'
gi=[
1 2 2 1 6
2 1 6 0 0
0 0 1 1 5
-1 -1 -1 -1 0
0 -1 0 -1 -1
-6 0 0 0 -7
1 1 1 1 1
-1 -1 -1 -1 -1
1 0 0 0 0
0 1 0 0 0
0 0 1 0 0
0 0 0 1 0
0 0 0 0 1
-1 0 0 0 0
0 -1 0 0 0
0 0 -1 0 0
0 0 0 -1 0
0 0 0 0 -1
];
alpha =size(gi,1)*1; % alpha-constraint
bb=[800,200,200,-48,-34,-104,400,-55, 99,99,99,99,99,0,0,0,0,0]';
b3=10;
kk= round(0.4*M);
delta=100;
pbest=zeros(2,M);
v=zeros(D+1,M);
globmin=10000000;
npr=0.1;
for i=1:1e7
    i
    s = rng;
    x=initialization(xmax,xmin);
    cnt=0;
    while delta>eps && cnt<10
       xbest=diff_evo(x);
       x=xbest;
       EDP(xbest);
       xx=findbestind(x);
       EDP(xx);
       if mod(g,200)==0
            pbest(1,1:M)=pbest(2,1:M);
            xbestind(1,1:M)=findbestind(x);
            xbestind(1,1:M)
            EDP(xbestind(1,1:M))
            P0=xbestind(1,1:M)';
            Pn=P0*npr*normrnd(mu,sigma,1,D);
            %Pn=P0*0.1*(rand(1,D)-0.5);
            Point=round(P0*ones(1,D)+Pn)';
            v(1,1:M)=P0;
            v(2:D+1,1:M)=Point;
            %fminsearch
            pbest(2,1:M)=discrete_nelder_mead(v)
            delta=abs(EDP(pbest(2,1:M))-EDP(pbest(1,1:M)));
            if getmu(pbest(2,1:M))<alpha
                delta=2*eps;
                cnt=cnt+1;
            end
       end
       g=g+1;
    end
    fmin=EDP(pbest(2,1:M))
    pbest(2,1:M)
    if((fmin<globmin) && (getmu(pbest(2,1:M))==alpha))
        globmin=fmin
        globmin
        globp=pbest(2,1:M);
        delta=100;
        cnt=0;
        while delta>eps && cnt<10
            P0=pbest(2,1:M)';
            Pn=P0*npr*normrnd(mu,sigma,1,D);
            %Pn=P0*0.1*(rand(1,D)-0.5);
            Point=round(P0*ones(1,D)+Pn)';
            v(1,1:M)=pbest(2,1:M);
            v(2:D+1,1:M)=Point;
            %fminsearch
            pbest(2,1:M)=discrete_nelder_mead(v)
            delta=abs(EDP(pbest(2,1:M))-EDP(P0'));
             if((EDP(pbest(2,1:M))<globmin) && (getmu(pbest(2,1:M))==alpha))
                globmin=EDP(pbest(2,1:M));
                globp=pbest(2,1:M);
             end
            if getmu(pbest(2,1:M))<alpha
                delta=2*eps;
                cnt=cnt+1;
            end
        end
    end
%    rng(s);
end
