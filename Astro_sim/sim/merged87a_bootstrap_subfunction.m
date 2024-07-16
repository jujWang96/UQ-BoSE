function [centx,centy,area,num_in,intensity,min_BIC,centx_seq,centy_seq,area_seq,num_in_B_seq,intensity_B_seq,num_in_O_seq,intensity_O_seq,BIC_seq]= merged87a_bootstrap_subfunction(X,target_curve,B,file_id,contour_merge,selected_merge)
%X =unique(X,'rows');
%keep original number of observations

    [centx,centy] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
    area = polyarea(target_curve(:,1),target_curve(:,2));
   
    [x_sample,y_sample] = sample_curve(target_curve(:,1),target_curve(:,2),300,centx, false);
    origin_coef = FD(x_sample,y_sample);
    minAIC = findminIC(X,origin_coef,'AIC',true,40);
    [invx,invy] = iFD(origin_coef,minAIC);
    num_in = sum(inpolygon(X(:,1),X(:,2),target_curve(:,1),target_curve(:,2)));
    intensity = num_in/area;
    %calculate the updated BIC value for two region segmentation
    min_BIC = calc_BIC_2seg(X,target_curve,10);
    seed = file_id;
    [bootstrap_X_set_aic,~] = sim_fix_data([0 1], [0 1],[invx,invy] ,sum(inpolygon(X(:,1),X(:,2),invx,invy))/length(X),length(X),B,seed);
    idx = 0;
    centx_seq=0;
    area_seq=0;
    num_in_B_seq=0;
    intensity_B_seq=0;
    num_in_O_seq=0;
    intensity_O_seq=0;
    BIC_seq=0;
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
            [centx_seq(idx),centy_seq(idx)] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
            area_seq(idx) = polyarea(target_curve(:,1),target_curve(:,2));
            num_in_B_seq(idx) = sum(inpolygon(X_boot(:,1),X_boot(:,2),target_curve(:,1),target_curve(:,2)));
            intensity_B_seq(idx) = num_in_B_seq(idx)/area_seq(idx);
            num_in_O_seq(idx) = sum(inpolygon(X(:,1),X(:,2),target_curve(:,1),target_curve(:,2)));
            intensity_O_seq(idx) = num_in_O_seq(idx)/area_seq(idx);
            BIC_seq(idx) = calc_BIC_2seg(X_boot,target_curve,10);
	end

    end
disp('number of valid segmentation is  ')
disp(idx)
end
    
