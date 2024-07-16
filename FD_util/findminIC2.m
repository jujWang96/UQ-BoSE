function minindex = findminIC2(X,f_coef,inner_bd_x,inner_bd_y, out_bd_x,out_bd_y,IC,pair,MM)

%find the number of fourier coefficients corresponds to minimum information
%criterion
%Input:
%   X: photon coordinates
%   IB_mat: Inner boundary of contour
%   OB: Outer boundary of contour
%   f_coef: K by d matrix. fourier descriptors of the contours
%   IC: type of criterion
%   MM: The maximum number of coefficients considered
%   pair: logic variable, if true then only consider pairs of coefficient. 

full_bd_x = inner_bd_x;
full_bd_y = inner_bd_y;
full_bd_x{end+1} = out_bd_x;
full_bd_y{end+1} = out_bd_y;
pgonRange = polyshape(full_bd_x,full_bd_y);
IB_mat=[];
if ~isempty(inner_bd_x)
    for d =1:length(inner_bd_x)
        curr_x = inner_bd_x{d};
        curr_y = inner_bd_y{d};
        IB_mat = [IB_mat;[curr_x,curr_y];[curr_x(1),curr_y(1)];NaN,NaN];
        
    end
    IB_mat(end,:)=[];

end

[~, ncol] = size(f_coef);
if ncol==1
    if strcmp(IC, 'BIC')
        if nargin == 8
            MM = 25;
        end
        Fbic = zeros(1,MM);
        for m = 3:2:MM
            Fbic(m) = FourierBIC(f_coef,m);
        end

        [minBIC,minindex] = min(Fbic);
    elseif strcmp(IC, 'AIC')
        if nargin == 8
            MM = 41;
        end
        Faic = zeros(1,MM);
        for m = 3:2:MM
        %for m = 13
            Faic(m) = FourierAIC(f_coef,m);
        end
        [minAIC,minindex] = min(Faic);

    end
else
    if strcmp(IC, 'BIC')
        if nargin == 8
            MM = 25;
        end
        minBIC = Inf;
        minindex = [0,0];
        for m = 3:2:MM
            for n = 3:2:MM
                currBIC = FourierBIC_two(f_coef,[m,n]);
                if currBIC < minBIC
                    minBIC = currBIC;
                    minindex = [m,n];
                end
                
            end
        end

    elseif strcmp(IC, 'AIC')
        if nargin == 8
            MM = 41;
        end
        minAIC = Inf;
        minindex = [0,0];
        for m = 3:2:MM
            for n = 3:2:MM
                currAIC = FourierAIC_two(f_coef,[m,n]);
                if currAIC < minAIC
                    minAIC = currAIC;
                    minindex = [m,n];
                end
                
            end
        end
    end
    
end
function bic = FourierBIC(f_coef, m)
    [xv,yv] = iFD(f_coef,m);
        xv(end+1) = xv(1);
    yv(end+1) = yv(1);
    pgon = polyshape(xv,yv);
    %drop the specific number of pairs if the recovered shape insect with
    %boundary or boundary cross
    if numboundaries(pgon)>1 || ~isempty(polyxpoly(xv,yv,out_bd_x,out_bd_y))...
        || ~(isempty(IB_mat)||isempty(polyxpoly(xv,yv,IB_mat(:,1),IB_mat(:,2)))) 
        bic=0;
        return 
    end
    
    if ~isempty(IB_mat)
        pgonIn = polyshape([IB_mat(:,1);NaN;xv],[IB_mat(:,2);NaN;yv]);
    else
        pgonIn = polyshape(xv,yv);
    end
    pgonOut = polyshape({out_bd_x;xv},{out_bd_y;yv});
    bic = -2*(sum(isinterior(pgonIn,X(:,1),X(:,2)))*log(sum(isinterior(pgonIn,X(:,1),X(:,2)))/area(pgonIn))...
        +sum(isinterior(pgonOut,X(:,1),X(:,2)))*log(sum(isinterior(pgonOut,X(:,1),X(:,2)))/area(pgonOut)))...
        +2*log(sum(isinterior(pgonRange,X(:,1),X(:,2))))*m;
    
end
function aic = FourierAIC(f_coef, m)
    [xv,yv] = iFD(f_coef,m);
    xv(end+1) = xv(1);
    yv(end+1) = yv(1);

    pgon = polyshape(xv,yv);
    %drop the specific number of pairs if the recovered shape insect with
    %boundary or boundary cross
    if numboundaries(pgon)>1 || ~isempty(polyxpoly(xv,yv,out_bd_x,out_bd_y))...
       || ~(isempty(IB_mat)||isempty(polyxpoly(xv,yv,IB_mat(:,1),IB_mat(:,2)))) 
        if (~isempty(polyxpoly(xv,yv,out_bd_x,out_bd_y)))
            disp('cross outer bound')
        end
        if ~(isempty(IB_mat)||isempty(polyxpoly(xv,yv,IB_mat(:,1),IB_mat(:,2)))) 
                        disp('cross inner bound')

        end
        aic=0;
        return 
    end
    
    if ~isempty(IB_mat)
        pgonIn = polyshape([IB_mat(:,1);NaN;xv],[IB_mat(:,2);NaN;yv]);
    else
        pgonIn = polyshape(xv,yv);
    end
    pgonOut = polyshape([out_bd_x;NaN;xv],[out_bd_y;NaN;yv]);
    plot(pgonIn)
    hold on
    plot(pgonOut)
    aic = -2*(sum(isinterior(pgonIn,X(:,1),X(:,2)))*log(sum(isinterior(pgonIn,X(:,1),X(:,2)))/area(pgonIn))...
        +sum(isinterior(pgonOut,X(:,1),X(:,2)))*log(sum(isinterior(pgonOut,X(:,1),X(:,2)))/area(pgonOut)))...
        +2*2*m;
    
    
end
function bic = FourierBIC_two(f_coef, num_coef)
    
[xv1,yv1] = iFD(f_coef(:,1),num_coef(:,1));
    xv1(end+1) = xv1(1);
    yv1(end+1) = yv12(1);
    pgon1 = polyshape(xv1,yv1);
    %drop the specific number of pairs if the recovered shape insect with
    %inner boundary or has a cross
    if numboundaries(pgon1)>1 || ~isempty(polyxpoly(xv1,yv1,out_bd_x,out_bd_y))...
        || ~(isempty(IB_mat)||isempty(polyxpoly(xv1,yv1,IB_mat(:,1),IB_mat(:,2)))) 
        bic=0;
        return 
    end
    if ~isempty(IB_mat) 
        in_bool = isinterior(pgon1,IB_mat(:,1),IB_mat(:,2));
        if sum(in_bool)>1
             pgonIn1 = polyshape([IB_mat(in_bool,1);NaN;xv1],[IB_mat(in_bool,2);NaN;yv1]);
        else
            pgonIn1 = pgon1;
        end
    else
        pgonIn1 = pgon1;
    end
    num1 = sum(isinterior(pgonIn1,X(:,1),X(:,2)));
    area1 = area(pgonIn1);
    [xv2,yv2] = iFD(f_coef(:,2),num_coef(:,2));
    pgon2 = polyshape(xv2,yv2);
    %drop the specific number of pairs if the recovered shape insect with
    %inner boundary or has a cross
   if numboundaries(pgon2)>1 || ~isempty(polyxpoly(xv2,yv2,out_bd_x,out_bd_y))...
        || ~(isempty(IB_mat)||isempty(polyxpoly(xv2,yv2,IB_mat(:,1),IB_mat(:,2)))) 
        bic=0;
        return 
    end
    if ~isempty(IB_mat) 
        in_bool = isinterior(pgon2,IB_mat(:,1),IB_mat(:,2));
        if sum(in_bool)>1
             pgonIn2 = polyshape([IB_mat(in_bool,1);NaN;xv2],[IB_mat(in_bool,2);NaN;yv2]);
        else
            pgonIn2 = pgon2;
        end
    else
        pgonIn2 = pgon2;
    end
     num2 = sum(isinterior(pgonIn2,X(:,1),X(:,2)));
    area2 = area(pgonIn2);
    
    pgonOut = polyshape([[out_bd_x;NaN;xv1;NaN;xv2],[out_bd_y;NaN;yv1;NaN;yv2]]);
    bic = -2*(num1*log(num1/area1)+num2*log(num2/area2)...
        +sum(isinterior(pgonOut,X(:,1),X(:,2)))*log(sum(isinterior(pgonOut,X(:,1),X(:,2)))/area(pgonOut)))...
        +2*log(sum(isinterior(pgonRange,X(:,1),X(:,2))))*sum(num_coef);
    
end
    function aic = FourierAIC_two(f_coef, num_coef)
    
   [xv1,yv1] = iFD(f_coef(:,1),num_coef(:,1));
    xv1(end+1) = xv1(1);
    yv1(end+1) = yv1(1);
    pgon1 = polyshape(xv1,yv1);
    %drop the specific number of pairs if the recovered shape insect with
    %inner boundary or has a cross
    if numboundaries(pgon1)>1 || ~isempty(polyxpoly(xv1,yv1,out_bd_x,out_bd_y))...
        || ~(isempty(IB_mat)||isempty(polyxpoly(xv1,yv1,IB_mat(:,1),IB_mat(:,2)))) 
        aic=0;
        return 
    end
    if ~isempty(IB_mat) 
        in_bool = isinterior(pgon1,IB_mat(:,1),IB_mat(:,2));
        if sum(in_bool)>1
             pgonIn1 = polyshape([IB_mat(in_bool,1);NaN;xv1],[IB_mat(in_bool,2);NaN;yv1]);
        else
            pgonIn1 = pgon1;
        end
    else
        pgonIn1 = pgon1;
    end
    num1 = sum(isinterior(pgonIn1,X(:,1),X(:,2)));
    area1 = area(pgonIn1);
    [xv2,yv2] = iFD(f_coef(:,2),num_coef(:,2));
    pgon2 = polyshape(xv2,yv2);
    %drop the specific number of pairs if the recovered shape insect with
    %inner boundary or has a cross
   if numboundaries(pgon2)>1 || ~isempty(polyxpoly(xv2,yv2,out_bd_x,out_bd_y))...
        || ~(isempty(IB_mat)||isempty(polyxpoly(xv2,yv2,IB_mat(:,1),IB_mat(:,2)))) 
        aic=0;
        return 
    end
    if ~isempty(IB_mat) 
        in_bool = isinterior(pgon2,IB_mat(:,1),IB_mat(:,2));
        if sum(in_bool)>1
             pgonIn2 = polyshape([IB_mat(in_bool,1);NaN;xv2],[IB_mat(in_bool,2);NaN;yv2]);
        else
            pgonIn2 = pgon2;
        end
    else
        pgonIn2 = pgon2;
    end
     num2 = sum(isinterior(pgonIn2,X(:,1),X(:,2)));
    area2 = area(pgonIn2);
    
    pgonOut = polyshape([[out_bd_x;NaN;xv1;NaN;xv2],[out_bd_y;NaN;yv1;NaN;yv2]]);

    aic = -2*(num1*log(num1/area1)+num2*log(num2/area2)...
        +sum(isinterior(pgonOut,X(:,1),X(:,2)))*log(sum(isinterior(pgonOut,X(:,1),X(:,2)))/area(pgonOut)))...
        +2*2*sum(num_coef);
    
end

end



