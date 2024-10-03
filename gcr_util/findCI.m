function CI = findCI(X,bootstrap_X_set, origin_coef, bootstrap_coef,Contour,num_coef,r,a,statistics)
%create the confidence band based on the 
%input:
%   bootstap_coef,num_coef: used to parametries the the boundaris of
%   bootstrapped data, bounded in [0 1]x[0 1]
%   r: number of grid in row and column 
%output:
%   CI: a rxr matrix represent the binaries grids. 1 if the grid is in side
%   the 1-a confidence band; 0 otherwise. 
%note: If the pixel's central subpixel is inside the modified polygon, the pixel is inside the region.
%the statistics could be profile loglikelihood or the difference of the
%shape 

n = length(num_coef);
nx = length(X);
%n = ceil(n*0.2);
hist = zeros(1,n);
holenum = 0;
if strcmp(statistics, 'logl')
    for l = 1:n
        %hist(l) =profilelogL(X,bootstrap_coef(:,l),num_coef(l))-log(nx)*num_coef(l);
        %hist(l) =profilelogL(X,bootstrap_coef(:,l),num_coef(l))-2*num_coef(l);
        
        hist(l) =profilelogL(X,bootstrap_coef(:,l),num_coef(l));
        [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
        pgon = polyshape(invx,invy);
        if numboundaries(pgon)>1
            holenum = holenum+1;
        end
        
    end
    disp('number of hole:')
    disp(holenum)
    %the empirical distribution is one-sided
    bd = prctile(hist,a*100);
    count = zeros(r,r);
    num=0;
    figure
    for l=1:n
        %if hist(l)>=prctile(hist,38) && hist(l)<=prctile(hist,38.5)
        if hist(l)>=bd
            num=num+1;
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            
            count = count+double(poly2mask(invx*r,invy*r,r,r));
        
        else
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            plot(invx,invy);
            axis equal
            hold on
        end  
    end
%use the loglikelihood ratio 
elseif strcmp(statistics, 'LR')
    l=0;
    for b = 1:length(bootstrap_X_set)
        if ~isempty(Contour{b})
            l=l+1;
            hist(l) =LRatio(bootstrap_X_set{b},origin_coef,bootstrap_coef(:,l),3,num_coef(l));

        end
        if l==n
            break
        end
    
    end
    %the empirical distribution is one-sided
    bd = prctile(hist,a*100);

    count = zeros(r,r);
    num=0;
    figure
    for l=1:n
        if hist(l)>=bd
            num=num+1;
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            plot(invx,invy);
            axis equal
            hold on
            count = count+double(poly2mask(invx*r,invy*r,r,r));
        end
    end
%use the distance between the two shape 
elseif strcmp(statistics, 'kld')

    for l = 1:n
        hist(l) =KLD(X,bootstrap_coef(:,l),origin_coef,num_coef(l),3);

    end
    %the empirical distribution is one-sided
    bd = prctile(hist,(1-a)*100);
    count = zeros(r,r);
    num=0;
    figure
    for l=1:n
        if hist(l)<=bd
            num=num+1;
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            
            count = count+double(poly2mask(invx*r,invy*r,r,r));
        
        else
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            plot(invx,invy);
            axis equal
            hold on
        end  
    end

%using size,center, shape seperately to decide 
else
    origin_coef_inv = make_inv(origin_coef,1);
    for l = 1:n
        %hist(l) =profilelogL(X,bootstrap_coef(:,l),num_coef(l))-log(nx)*num_coef(l);
        %hist(l) =profilelogL(X,bootstrap_coef(:,l),num_coef(l))-2*num_coef(l);
        histct(l) =dist_complex(origin_coef,bootstrap_coef(:,l),1,false);
        histsz(l) = abs(abs(origin_coef(2))-abs(bootstrap_coef(2,l)));
        [coef_a,coef_b] = make_inv(bootstrap_coef(:,l),floor((num_coef-1)/2));
        histsp(l)= min(dist_complex(coef_a,origin_coef_inv),dist_complex(coef_b,origin_coef_inv));

        
    end
    %the empirical distribution is one-sided
    bdct = prctile(histct,(1-a/3)*100);
    bdsz = prctile(histsz,(1-a/3)*100);
    bdsp = prctile(histsp,(1-a/3)*100);

    count = zeros(r,r);
    num=0;
    figure
    for l=1:n
        
        if histct(l)<bdct && histsz(l)<=bdsz && histsp(l)<=bdsp
            num=num+1;
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            
            count = count+double(poly2mask(invx*r,invy*r,r,r));
        
        else
            [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
            plot(invx,invy);
            axis equal
            hold on
        end  
    end
end
disp(num)
CI = zeros(r,r);
CI(count>0&count<num)=1;
CI= flip(CI);

