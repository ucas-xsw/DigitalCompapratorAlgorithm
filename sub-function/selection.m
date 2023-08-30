function xbest=selection(x,u)
xbest=x;
global NP
    for i=1:NP
        select=alphaleq(x(i,:),u(i,:));
        if select==2
            xbest(i,:)=u(i,:);
        end
    end

end