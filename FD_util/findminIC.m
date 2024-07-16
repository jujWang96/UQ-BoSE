function minindex = findminIC(X,f_coef,IC,pair,MM)
%find the number of fourier coefficients corresponds to minimum information
%criterion
%Input:
%MM: The maximum number of coefficients considered
%pair: logic variable, if true then only consider pairs of coefficient. 

[~,d] = size(f_coef);

if pair == false
    if strcmp(IC, 'AIC')
        if nargin == 4
            MM = 20;
        end
        if d==1
            Faic = zeros(1,MM);
            for m = 1:MM
                Faic(m) = FourierAIC(X,f_coef,m);
            end
            [minAIC,minindex] = min(Faic);
        else
            Faic = zeros(MM,MM);
            for i = 1:MM
                for j = 1:MM
                    Faic(i,j) = FourierAIC(X,f_coef,[i,j]);
                end
            end
            [min_val,idx]=min(Faic(:));
            [row,col]=ind2sub(size(Faic),idx);
            minindex = [row,col];
        end

    elseif strcmp(IC,'BIC') 
        if nargin == 4
            MM = 10;
        end
        Fbic = zeros(1,MM);
        for m = 1:MM
            Fbic(m) = FourierBIC(X,f_coef,m);
        end
        [minBIC,minindex] = min(Fbic);
    elseif strcmp(IC,'AICc') 
        if nargin == 4
            MM = 15;
        end
        Faicc = zeros(1,MM);
        for m = 1:MM
            Faicc(m) = FourierAICc(X,f_coef,m);
        end
        [minAICc,minindex] = min(Faicc);
    elseif strcmp(IC,'MDL')
        if nargin == 4
            MM = 20;
        end
        Fmdl = zeros(1,MM);
        for m = 1:MM
            Fmdl(m) = FourierMDL(X,f_coef,m);
        end
        [minMDL,minindex] = min(Fmdl);
    else
        if nargin == 4
            MM = 15;
        end
        Fhqc = zeros(1,MM);
        for m = 1:MM
            Fhqc(m) = FourierHQC(X,f_coef,m);
        end
        [minHQC,minindex] = min(Fhqc);
    end
else
    if strcmp(IC, 'AIC')
        if nargin == 4
            MM = 20;
        end
         if d==1
            Faic = zeros(1,MM);
            for m = 1:2:MM
                Faic(m) = FourierAIC(X,f_coef,m);
            end
            %Faic = nonzeros(Faic);
            [minAIC,minindex] = min(Faic);
            disp(Faic)
            %minindex = minindex*2-1;
        else
            Faic = zeros(MM,MM);
            for i = 1:2:MM
                for j = 1:2:MM
                    Faic(i,j) = FourierAIC(X,f_coef,[i,j]);
                end
            end
            [min_val,idx]=min(Faic(:));
            
            [row,col]=ind2sub(size(Faic),idx);
            minindex = [row,col];
        end
    elseif strcmp(IC,'BIC') 
        if nargin == 4
            MM = 10;
        end
        Fbic = zeros(1,MM);
        for m = 1:2:MM
            Fbic(m) = FourierBIC(X,f_coef,m);
        end
        [minBIC,minindex] = min(Fbic);

    elseif strcmp(IC,'AICc') 
        if nargin == 4
            MM = 15;
        end
        Faicc = zeros(1,MM);
        for m = 1:2:MM
            Faicc(m) = FourierAICc(X,f_coef,m);
        end
        Faicc = nonzeros(Faicc);
        [minAICc,minindex] = min(Faicc);
     elseif strcmp(IC,'MDL')
        if nargin == 4
            MM = 20;
        end
        Fmdl = zeros(1,MM);
        for m = 1:MM
            Fmdl(m) = FourierMDL(X,f_coef,m);
        end
        Fmdl = nonzeros(Fmdl);
        [minMDL,minindex] = min(Fmdl);
    else
        if nargin == 4
            MM = 15;
        end
        Fhqc = zeros(1,MM);
        for m = 1:1:MM
            Fhqc(m) = FourierHQC(X,f_coef,m);
        end
        Fhqc = nonzeros(Fhqc);
        [minHQC,minindex] = min(Fhqc);
        minindex = minindex*2-1;
    end
end



