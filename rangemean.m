function  S=rangemean(xs,x,y,twin,plotflag)
% [mx,my,myse,mxse,mxvar,myvar,mykurt]=rangemean(xs,x,y,twin,plotflag)
%   Finds the mean (and standard errors) of the b points spanning the
%   sample points xs from the (x,y) data. At the edges of the data only b/2
%   points may be averaged - and optionally plots the data



for i=1:length(xs)
    k=find(x>=xs(i)-twin & x<=xs(i)+twin);
    
    if isempty(k)
        xs(i)
        pause
    end
    
    
    %ensure window does not extend past data
    
    mx(i)=mean(x(k));
    my(i)=mean(y(k));
    mxse(i)=stderr(x(k));
    myse(i)=stderr(y(k));
    mxvar(i)=var(x(k));
    myvar(i)=var(y(k));
    
    mykurt(i)=kurtosis(y(k));
    
end

if nargin ==5 & plotflag
    plot(mx,my)
    hold on
    plot(mx,my+myse)
    plot(mx,my-myse)
    shg
end

S=v2struct(mx,my,myse,mxse,mxvar,myvar,mykurt);

