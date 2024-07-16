function [centx,centy,area,num_in,intensity,min_BIC,centx_seq,centy_seq,area_seq,num_in_B_seq,intensity_B_seq,num_in_O_seq,intensity_O_seq,BIC_seq,BIC_O_seq,centx_fd,centy_fd,area_fd,num_in_fd,intensity_fd,min_BIC_fd,centx_seq_fd,centy_seq_fd,area_seq_fd,num_in_B_seq_fd,intensity_B_seq_fd,num_in_O_seq_fd,intensity_O_seq_fd,BIC_seq_fd,BIC_O_seq_fd]= merged87a_bootstrap_smooth_subfunction(X,target_curve,B,file_id,tline)
%X =unique(X,'rows');
%keep original number of observations
    %calculate statistics of raw segmentation result without model selection 
    [centx,centy] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
    num_in = sum(inpolygon(X(:,1),X(:,2),target_curve(:,1),target_curve(:,2)));
    area = polyarea(target_curve(:,1),target_curve(:,2));
    intensity = num_in/area;
    min_BIC = calc_BIC_2seg(X,target_curve,10);
    %calculate statistics of smoothed curve with model selection  
    [x_sample,y_sample] = sample_curve(target_curve(:,1),target_curve(:,2),300,centx, false);
    origin_coef = FD(x_sample,y_sample);
    centx_fd = real(origin_coef(1));
    centy_fd = imag(origin_coef(1));
    min_modelselect = findminIC(X,origin_coef,'AIC',true,40);
    [invx,invy] = iFD(origin_coef,min_modelselect);
    area_fd = polyarea(invx,invy);
    in_seg = inpolygon(X(:,1),X(:,2),invx,invy);
    num_in_fd = sum(in_seg);
    intensity_fd = num_in_fd/area_fd;
    %calculate the updated BIC value for two region segmentation
    min_BIC_fd = calc_BIC_2seg(X,[invx,invy],10);
    seed = file_id;
    [bootstrap_X_set_aic,N_in_out] = sim_fix_data([0 1], [0 1],[invx,invy] ,sum(in_seg)/length(X),length(X),B,seed);
    N_in_out
    idx = 0;
    centx_seq=0;
    centy_seq=0;
    area_seq=0;
    num_in_B_seq=0;
    intensity_B_seq=0;
    num_in_O_seq=0;
    intensity_O_seq=0;
    BIC_seq=0;
    BIC_O_seq = 0;
    centx_seq_fd=0;
    centy_seq_fd=0;
    area_seq_fd=0;
    num_in_B_seq_fd=0;
    intensity_B_seq_fd=0;
    num_in_O_seq_fd=0;
    intensity_O_seq_fd=0;
    BIC_seq_fd=0;
    BIC_O_seq_fd = 0;
    for b =1:B
        X_boot = bootstrap_X_set_aic{b};
        [boot_contour,boot_region_intensity,boot_selected,boot_min_BIC] = gSRG_random3(double(unique(X_boot,'rows')),false,false,'k',[],true,true,2,30,10000,10);
	%[boot_contour,boot_selected] = segment_graph_contour(X_boot,false,false, 'k',1,true,true,10)
	objId = setdiff(boot_contour.regionId,boot_contour.background);
        if length(objId)>1 || length(objId)<1
            continue
        else
            idx=idx+1;
            clear target_curve
        
            [target_curve(:,1),target_curve(:,2)] = get_curve(boot_contour.contourV{objId},false);    
            [centx,centy] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
             %calculate bootstrap statistics of raw segmentation without model selection
            
            centx_seq(idx) = centx;
            centy_seq(idx) = centy;
            area_seq(idx) = polyarea(target_curve(:,1),target_curve(:,2));
            num_in_B_seq(idx) = sum(inpolygon(X_boot(:,1),X_boot(:,2),target_curve(:,1),target_curve(:,2)));
            intensity_B_seq(idx) = num_in_B_seq(idx)/area_seq(idx);
            num_in_O_seq(idx) = sum(inpolygon(X(:,1),X(:,2),target_curve(:,1),target_curve(:,2)));
            intensity_O_seq(idx) = num_in_O_seq(idx)/area_seq(idx);
            BIC_O_seq(idx) = calc_BIC_2seg(X,target_curve,10);

            BIC_seq(idx) = calc_BIC_2seg(X_boot,target_curve,10);
            %calculate bootstrap statistics after model selection
            [x_sample,y_sample] = sample_curve(target_curve(:,1),target_curve(:,2),300,centx, false);
            origin_coef = FD(x_sample,y_sample);
            min_modelselect = findminIC(X,origin_coef,'AIC',true,40);
            [invx,invy] = iFD(origin_coef,min_modelselect);
            centx_seq_fd(idx) = real(origin_coef(1));
            centy_seq_fd(idx) = imag(origin_coef(1));
            area_seq_fd(idx) = polyarea(invx,invy);
            num_in_B_seq_fd(idx) = sum(inpolygon(X_boot(:,1),X_boot(:,2),invx,invy));
            intensity_B_seq_fd(idx) = num_in_B_seq_fd(idx)/area_seq_fd(idx);
            num_in_O_seq_fd(idx) = sum(inpolygon(X(:,1),X(:,2),invx,invy));
            intensity_O_seq_fd(idx) = num_in_O_seq_fd(idx)/area_seq_fd(idx);
            BIC_O_seq_fd(idx) = calc_BIC_2seg(X,[invx,invy],10);

            BIC_seq_fd(idx) = calc_BIC_2seg(X_boot,[invx,invy],10);
            boot_contours{idx} = boot_contour;
            boot_region_intensitys{idx} = boot_region_intensity;
            boot_selecteds{idx} = boot_selected;
        end

    end
save(strcat('1987a_bootstrap_FD/',tline,'bootstrap_segmentation.mat'),'boot_contours','boot_region_intensitys','boot_selecteds','bootstrap_X_set_aic')
disp('number of valid segmentation is  ')
disp(idx)
end
    
