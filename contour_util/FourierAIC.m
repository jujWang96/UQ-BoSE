function  aic = FourierAIC(X,coef,M)
%
%X: input dataset
%coef: K by d matrix. Fourier coefficients of d segmentation contours
%assuming d=2 at most at this time point 
%M: number of Fourier coefficients used to recover the shape
%N: numebr of vextices in recovered shape

d = length(M);
N = length(X);
cx = X(:,1);
cy = X(:,2);

if d==1
    [xv,yv] = iFD(coef,M);
    pgon = polyshape(xv,yv);
% 
    if numboundaries(pgon)>1
        aic=0;
        return 
    end
    %  
    A = polyarea(xv,yv);
    in = inpolygon(cx,cy,xv,yv);
    out = ~in;
    aic = -2*(sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A)))+ 2*2*M;
else
    [xv1,yv1] = iFD(coef(:,1),M(1));
    [xv2,yv2] = iFD(coef(:,2),M(2));

    itsct = intersect(polyshape(xv1,yv1),polyshape(xv2,yv2));
    xvits= itsct.Vertices(:,1);
    yvits= itsct.Vertices(:,2);
    num1 = sum(inpolygon(cx,cy,xv1,yv1));
    num2 = sum(inpolygon(cx,cy,xv2,yv2));
    numits = sum(inpolygon(cx,cy,xvits,yvits));
    A1 = polyarea(xv1,yv1);
    A2 = polyarea(xv2,yv2);
    Aits = polyarea(xvits,yvits);
    aic = -2*((num1-numits)*log((num1-numits)/(A1-Aits))+...
        (num2-numits)*log((num2-numits)/(A1-Aits))+...
        numits*log(numits)/Aits+...
        (N-num1-num2+numits)*log((N-num1-num2+numits)/(1-A1-A2+Aits)))+ 2*2*(M(1)+M(2));
end





