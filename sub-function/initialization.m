function x=initialization(xmax,xmin)
global NP M
x= ones(NP,M)*xmin+rand(NP,M)*(xmax-xmin);
end