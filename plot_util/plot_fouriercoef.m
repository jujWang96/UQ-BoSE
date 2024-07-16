addpath(genpath('/Users/joycewang/Desktop/segmentation/G-SRG'))
addpath(genpath('/Users/joycewang/Desktop/segmentation/Astro_sim'))


%load boot_circle_1000.mat
load check_X_018_si.mat
load boot_circle_1000.mat
F_coef = boot_F_coef;
%load X_circle.mat
L = 200;
for l = 1:L
    scatter(real(F_coef(1,l)),imag(F_coef(1,l)),'.','r');
    hold on;
    scatter(real(F_coef(2,l)),imag(F_coef(2,l)),'.','b');
    hold on;
    scatter(real(F_coef(300,l)),imag(F_coef(300,l)),'.','m');
    hold on;
    %scatter(real(F_coef(3,k)),imag(F_coef(3,k)),'.','g');
    %hold on;
    %scatter(real(F_coef(299,k)),imag(F_coef(299,k)),[],[0.61 0.51 0.74],'.');
    %hold on;
    %scatter(real(F_coef(4,k)),imag(F_coef(4,k)),'.','c');
    %hold on;
    %scatter(real(F_coef(298,k)),imag(F_coef(298,k)),'.','y');
    %hold on;

end
 scatter(real(origin_coef(1)),imag(origin_coef(1)),'k','filled');
 text(real(origin_coef(1))+0.01,imag(origin_coef(1))+0.01,num2str(0));
    hold on;
    scatter(real(origin_coef(2)),imag(origin_coef(2)),'k','filled');
     text(real(origin_coef(2))+0.01,imag(origin_coef(2))+0.01,num2str(1));

    hold on;
    scatter(real(origin_coef(300)),imag(origin_coef(300)),'k','filled');
     text(real(origin_coef(300))+0.01,imag(origin_coef(300))+0.01,num2str(1));

    hold on;
    scatter(real(origin_coef(3)),imag(origin_coef(3)),'k','filled');
     text(real(origin_coef(3))+0.01,imag(origin_coef(3))+0.01,num2str(2));

    hold on;
    scatter(real(origin_coef(299)),imag(origin_coef(299)),'k','filled');
     text(real(origin_coef(299))+0.01,imag(origin_coef(299))+0.01,num2str(2));

    hold on;
    scatter(real(origin_coef(4)),imag(origin_coef(4)),'k','filled');
     text(real(origin_coef(4))+0.01,imag(origin_coef(4))+0.01,num2str(3));

    hold on;
    scatter(real(origin_coef(298)),imag(origin_coef(298)),'k','filled');
     text(real(origin_coef(298))+0.01,imag(origin_coef(298))+0.01,num2str(3));

    hold on;