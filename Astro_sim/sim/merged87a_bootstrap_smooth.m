fid = fopen('87a_mergeROIs_vlk_subset1.txt');
%outputfile = fopen('87a_vlk_output_2seg_subset1.txt','w');
B=5;
outputfile = fopen('experiment.txt','w');
%while ischar(tline)
for i = 1:1
    tline = fgetl(fid);
    clear target_curve
    try
        tline ='01044	41 90 92 107 108';
        read_line= str2num(tline);
        file_id = read_line(1)
        merge_seg = read_line(2:end);
        %disp(tline)
        fprintf(outputfile,'%s\n',num2str(read_line));
        load(strcat('1987a/acisf',num2str(file_id,'%05.f'),'_repro_evt2_scaled.mat'))
        load(strcat('1987a_randmerge/acisf',num2str(file_id,'%05.f'),'_rand3_1010_newseeds.mat'))
        [contour_merge,selected_merge] = merge_segments(contour,selected,merge_seg);
        contourV = contour_merge.contourV;
        [target_curve(:,1),target_curve(:,2)] = get_curve(contourV{merge_seg(1)},false);
        [centx,centy,area,num_in,intensity,min_BIC,centx_seq,centy_seq,area_seq,num_in_B_seq,intensity_B_seq,num_in_O_seq,intensity_O_seq,BIC_seq,BIC_O_seq,centx_fd,centy_fd,area_fd,num_in_fd,intensity_fd,min_BIC_fd,centx_seq_fd,centy_seq_fd,area_seq_fd,num_in_B_seq_fd,intensity_B_seq_fd,num_in_O_seq_fd,intensity_O_seq_fd,BIC_seq_fd,BIC_O_seq_fd] = merged87a_bootstrap_smooth_subfunction(X,target_curve,B,file_id,tline);
        fprintf(outputfile,'%s %f %f %f %f %f\n','centx',centx,[centx-std(centx_seq),centx+std(centx_seq),quantile(centx_seq,0.16),quantile(centx_seq,0.84)]);
        fprintf(outputfile,'%s %f %f %f %f %f\n','centy',centy,[centy-std(centy_seq),centy+std(centy_seq),quantile(centy_seq,0.16),quantile(centy_seq,0.84)]);
        fprintf(outputfile,'%s %f %f %f %f %f\n','area',area,[area-std(area_seq),area+std(area_seq),quantile(area_seq,0.16),quantile(area_seq,0.84)]);
        fprintf(outputfile,'%s %f %f %f %f %f\n','photon_in',num_in,[num_in-std(num_in_B_seq),num_in+std(num_in_B_seq),quantile(num_in_B_seq,0.16),quantile(num_in_B_seq,0.84)]);
	  fprintf(outputfile,'%s %f %f %f %f %f\n','photon_in_OBS',num_in,[num_in-std(num_in_O_seq),num_in+std(num_in_O_seq),quantile(num_in_O_seq,0.16),quantile(num_in_O_seq,0.84)]);
        fprintf(outputfile,'%s %f %f %f %f %f\n','brightness',intensity,[intensity-std(intensity_B_seq),intensity+std(intensity_B_seq),quantile(intensity_B_seq,0.16),quantile(intensity_B_seq,0.84)]);
        fprintf(outputfile,'%s %f %f %f %f %f\n','brightness_obs',intensity,[intensity-std(intensity_O_seq),intensity+std(intensity_O_seq),quantile(intensity_O_seq,0.16),quantile(intensity_O_seq,0.84)]);
	fprintf(outputfile,'%s %f %f %f %f %f\n','minBIC',min_BIC,[min_BIC-std(BIC_seq),min_BIC+std(BIC_seq),quantile(BIC_seq,0.16),quantile(BIC_seq,0.84)]);
    	fprintf(outputfile,'%s %f %f %f %f %f\n','minBIC',min_BIC,[min_BIC-std(BIC_O_seq),min_BIC+std(BIC_O_seq),quantile(BIC_O_seq,0.16),quantile(BIC_O_seq,0.84)]);
    
        
        fprintf(outputfile2,'%s %f %f %f %f %f\n','centx',centx_fd,[centx_fd-std(centx_seq_fd),centx_fd+std(centx_seq_fd),quantile(centx_seq_fd,0.16),quantile(centx_seq_fd,0.84)]);
        fprintf(outputfile2,'%s %f %f %f %f %f\n','centy',centy_fd,[centy_fd-std(centy_seq_fd),centy_fd+std(centy_seq_fd),quantile(centy_seq_fd,0.16),quantile(centy_seq_fd,0.84)]);
        fprintf(outputfile2,'%s %f %f %f %f %f\n','area',area,[area_fd-std(area_seq_fd),area+std(area_seq_fd),quantile(area_seq_fd,0.16),quantile(area_seq_fd,0.84)]);
        fprintf(outputfile2,'%s %f %f %f %f %f\n','photon_in',num_in_fd,[num_in_fd-std(num_in_B_seq_fd),num_in+std(num_in_B_seq_fd),quantile(num_in_B_seq_fd,0.16),quantile(num_in_B_seq_fd,0.84)]);
	  fprintf(outputfile2,'%s %f %f %f %f %f\n','photon_in_OBS',num_in_fd,[num_in_fd-std(num_in_O_seq_fd),num_in+std(num_in_O_seq_fd),quantile(num_in_O_seq_fd,0.16),quantile(num_in_O_seq_fd,0.84)]);
        fprintf(outputfile2,'%s %f %f %f %f %f\n','brightness',intensity_fd,[intensity_fd-std(intensity_B_seq_fd),intensity_fd+std(intensity_B_seq_fd),quantile(intensity_B_seq_fd,0.16),qunsantile(inteity_B_seq_fd,0.84)]);
        fprintf(outputfile2,'%s %f %f %f %f %f\n','brightness_obs',intensity_fd,[intensity_fd-std(intensity_O_seq_fd),intensity_fd+std(intensity_O_seq_fd),quantile(intensity_O_seq_fd,0.16),quantile(intensity_O_seq_fd,0.84)]);
	fprintf(outputfile2,'%s %f %f %f %f %f\n','minBIC',min_BIC_fd,[min_BIC_fd-std(BIC_seq_fd),min_BIC_fd+std(BIC_seq_fd),quantile(BIC_seq_fd,0.16),quantile(BIC_seq_fd,0.84)]);
    	fprintf(outputfile2,'%s %f %f %f %f %f\n','minBIC',min_BIC_fd,[min_BIC_fd-std(BIC_O_seq_fd),min_BIC_fd+std(BIC_O_seq_fd),quantile(BIC_O_seq_fd,0.16),quantile(BIC_O_seq_fd,0.84)]);
    
    fprintf(outputfile,'%s %f\n','photon_total',length(X));
        save(strcat('1987a_bootstrap_FD/',tline,'bootstrap_statistics.mat'),'centx','centy','area','num_in','intensity','min_BIC','centx_seq','centy_seq','area_seq','num_in_B_seq','intensity_B_seq','num_in_O_seq','intensity_O_seq','BIC_seq','centx_fd','centy_fd','area_fd','num_in_fd','intensity_fd','min_BIC_fd','centx_seq_fd','centy_seq_fd','area_seq_fd','num_in_B_seq_fd','intensity_B_seq_fd','num_in_O_seq_fd','intensity_O_seq_fd','BIC_seq_fd','BIC_O_seq_fd');   

    catch
        disp('exception: merged region cannot form a complete region')
    end
end
fclose(outputfile);

fclose(fid);
