clear 
cent = [0.5:0.01:0.56];
radi = [0.17:0.01:0.24];
powerlist = zeros(1,numel(radi));
powerlist_3 = zeros(1,numel(cent));
powerlist_11 = zeros(1,numel(cent));
powerlist_all = zeros(1,numel(cent));



%load check_result_020_first_90.mat
powerlist_area(2) = power;

%load check_result_010_si_first_90.mat
powerlist_area_si(1)=power;
LRpower3 = LRpower3(1:7);

powerlist_area_si(6) = powerlist_area(6); 
%datacursormode(h,'on');
plot(cent,LRpower1,'-o');
hold on

plot(cent,LRpower2,'-o');
hold on
plot(cent,LRpower3,'-o');


xlabel("center")
ylabel("power")
legend("no selection","BIC","AIC",'Location','southeast')
title("significant level=0.1, circle of different center")
for k = 1:numel(cent)

    text(cent(k),Bonfpower_ct(k),[num2str(round(Bonfpower_ct(k),3))])


end

save temp.mat powerlist_0 powerlist_3 powerlist_11
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

subplot(1,2,2);

load X_circle.mat

%load sim2_circle_1000.mat
load check_X_020.mat

i=0;
B=200;
for b = 1:B
    
    if ~isempty(Contour{b})
        i=i+1;
        if i == 31
            break;
        end
        
        contour = Contour{b};
        x1 = contour(:,1);
        x2 = contour(:,2);
        y1 = contour(:,3);
        y2 = contour(:,4);
        [x,y] = get_curve(x1,x2,y1,y2,true);
        [centx,centy] = centroid(polyshape(x,y));
        [x_sample,y_sample] = sample_curve(x,y,K,centx, false);
        
        %F_coef(:,i) = make_scale_inv(FD(x_sample,y_sample));
        coef = FD(x_sample,y_sample);
        coef = [0,coef(2),zeros(1,K-2)];
        a = ifft(coef*length(coef));
        
        invx = real(a);
        invy = imag(a);
        invx = [invx,invx(1)];
        invy = [invy,invy(1)];
        axis square;
        plot(invx,invy,'r');
        hold on;
  
    end
    
end

load check_X_022_si.mat
B=200;
i=0;
for b = 1:B
    if ~isempty(Contour{b})
        i=i+1;
        if i == 21
            break;
        end
        
        contour = Contour{b};
        x1 = contour(:,1);
        x2 = contour(:,2);
        y1 = contour(:,3);
        y2 = contour(:,4);
        [x,y] = get_curve(x1,x2,y1,y2,true);
        [centx,centy] = centroid(polyshape(x,y));
        [x_sample,y_sample] = sample_curve(x,y,K,centx, false);
        
        %F_coef(:,i) = make_scale_inv(FD(x_sample,y_sample));
        coef = FD(x_sample,y_sample);
        coef = [0,coef(2),zeros(1,K-2)];
        a = ifft(coef*length(coef));
        
        invx = real(a);
        invy = imag(a);
        invx = [invx,invx(1)];
        invy = [invy,invy(1)];
        axis square;
        plot(invx,invy,'b');
        hold on;
  
    end
    
end

load check_X_018_si.mat
B=200;
i=0;
for b = 1:B
    if ~isempty(Contour{b})
        i=i+1;
        if i == 21
            break;
        end
        
        contour = Contour{b};
        x1 = contour(:,1);
        x2 = contour(:,2);
        y1 = contour(:,3);
        y2 = contour(:,4);
        [x,y] = get_curve(x1,x2,y1,y2,true);
        [centx,centy] = centroid(polyshape(x,y));
        [x_sample,y_sample] = sample_curve(x,y,K,centx, false);
        
        %F_coef(:,i) = make_scale_inv(FD(x_sample,y_sample));
        coef = FD(x_sample,y_sample);
        coef = [0,coef(2),zeros(1,K-2)];
        a = ifft(coef*length(coef));
        
        invx = real(a);
        invy = imag(a);
        invx = [invx,invx(1)];
        invy = [invy,invy(1)];
        axis square;
        plot(invx,invy,'g');
        hold on;
  
    end
    
end

coef = [0,origin_coef(2),zeros(1,K-2)];
a = ifft(coef*length(coef));

invx = real(a);
invy = imag(a);
invx = [invx,invx(1)];
invy = [invy,invy(1)];
axis square;
plot(invx,invy,'k','LineWidth',2);
title(".18-.20-.22 bootstrap")


th = 0:2*pi/300:2*pi;

plot(cos(th)*0.2,sin(th)*0.2,'k','LineWidth',2);
hold on;
title(".18-.20-.22 groundtruth")

load sim2_dist_1000_first.mat
subplot(1,2,1);
histogram(distC);
title('histogram of bootstrap distance')
load groundtruth400_first.mat
subplot(1,2,2);
histogram(distC);
title('histogram of true distance')

K=300;
load check_X_020.mat 
invx = zeros(1,100);
for i =1:100
    coef = F_coef(:,i);
    invx(i) = abs(coef(2));

end

histogram(invx,0.16:0.005:0.23,'facealpha',0.3,'facecolor','b');
hold on

load boot_circle_1000.mat 
invx = zeros(1,100);
for i =1:100
    coef = F_coef(:,i);
    invx(i) = abs(coef(2));

end

histogram(invx,0.16:0.005:0.23,'facealpha',0.3,'facecolor','r');
hold on

load old_sim2/sim2_circle_1000.mat
invx = zeros(1,100);
for i =1:100
    coef = F_coef(:,i);
    invx(i) = abs(coef(2));

end

histogram(invx,0.16:0.005:0.23,'facealpha',0.3,'facecolor','g');
hold on


legend('ground truth','new bootstraping method','old boostrapping method')
title(' recovered circle radius')




subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);
load X_circle.mat
subplot(2,3,1)
scatter(X(:,1),X(:,2),'.');
axis([0 1 0 1])
axis square
hold on;
[invx,invy] = iFD(origin_coef);
plot(invx,invy)
title('circle 0.2')

load check_X_020.mat
subplot(2,3,2)
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),'.');
axis([0 1 0 1])
axis square
title('circle 0.2')
subplot(2,3,3)

load boot_circle_1000.mat
X_boot = bootstrap_X_set{1};
scatter(X_boot(:,1),X_boot(:,2),'.');
axis([0 1 0 1])
axis square
title('bootstrapped circle 0.2')


%load check_X_ellipse_032.mat
load check_X_018.mat
X_check = check_X_set{1};
subplot(2,3,4)
scatter(X_check(:,1),X_check(:,2),'.');
axis([0 1 0 1])
axis square
title('circle 0.18')


subplot(2,3,5)

load check_X_pois_022.mat
X_check = check_X_set{3};
scatter(X_check(:,1),X_check(:,2),'.');
axis([0 1 0 1])
axis square
title('circle 0.22')

load check_X_ellipse_018.mat
X_check = check_X_set{1};
subplot(2,3,6)
scatter(X_check(:,1),X_check(:,2),'.');
axis([0 1 0 1])
axis square
title('ellipse 0.18')


load check_X_ellipse_016_piby8.mat
X_check = check_X_set{3};
Fcof = F_coef(:,2);

[invx, invy] = iFD(Fcof,1);
subplot(1,2,1)
segment_graph_contour(X_check,true,true,'r',2,true,2);
subplot(1,2,2)
scatter(X_check(:,1),X_check(:,2),'.');
axis([0,1,0,1])
axis equal;
hold on;
plot(invx,invy);

for b = 1:100
    scatter
end
