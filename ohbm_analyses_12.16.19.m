% load the necessary data

%----- USING ROIS separately -----

% MAPS DATA
% Fit general linearl model
for i = 1:length(ROI_list)
    % Make a temporary table
    temp_T = [maps_thickness(:,ROI_list{i}) maps_info(:,{'age','sex','ari'})];
    lm_age_sex_ari_maps{i} = fitlm(temp_T,[ROI_list{i,:},'~age+sex+ari']);
    uncorrected_pvals(i) = lm_age_sex_ari_maps{i}.Coefficients{4,4};
end

% FDR correction
[h crit_p adj_ci adj_pvals_ari_maps] = fdr_bh(uncorrected_pvals);

% Display significant
ROI_list(h)

% TEMPER LOSS
% Fit general linearl model
for i = 1:length(ROI_list)
    % Make a temporary table
    temp_T = [maps_thickness(:,ROI_list{i}) maps_info(:,{'age','sex','temper_loss'})];
    lm_age_sex_temper_maps{i} = fitlm(temp_T,[ROI_list{i,:},'~age+sex+temper_loss']);
    uncorrected_pvals(i) = lm_age_sex_temper_maps{i}.Coefficients{4,4};
end

% FDR correction
[h crit_p adj_ci adj_pvals_temper_maps] = fdr_bh(uncorrected_pvals);

% Display significant
ROI_list(h)

% BDS DATA
% Fit general linear model
for i = 1:length(ROI_list)
    % Make a temporary table
    temp_T = [bds_thickness(:,ROI_list{i}) bds_info(:,{'MRIage','sex','ARIcat'})];
    lm_age_sex_ari_bds{i} = fitlm(temp_T,[ROI_list{i,:},'~MRIage+sex+ARIcat']);
    uncorrected_pvals(i) = lm_age_sex_ari_bds{i}.Coefficients{4,4};
end

% FDR correction
[h crit_p adj_ci adj_pvals_ari_bds] = fdr_bh(uncorrected_pvals);

% Display significant
ROI_list(h)

% COMBINED DATASETS
for i = 1:length(ROI_list)
    % Make a temporary table
    temp_T_bds = [bds_thickness(:,ROI_list{i}) bds_info(:,{'MRIage','sex','ARIcat'})];
    temp_T_bds.Properties.VariableNames{'MRIage'} = 'age';
    temp_T_bds.Properties.VariableNames{'ARIcat'} = 'ari';
    temp_T_maps = [maps_thickness(:,ROI_list{i}) maps_info(:,{'age','sex','ari'})];
    
    temp_T_total = [temp_T_maps;temp_T_bds];
    site = categorical([ones(size(temp_T_maps,1),1);zeros(size(temp_T_bds,1),1)]);
    temp_T_total = addvars(temp_T_total,site);

    lm_age_sex_ari_site_total{i} = fitlm(temp_T_total,[ROI_list{i,:},'~age+sex+site+ari']);
    uncorrected_pvals(i) = lm_age_sex_ari_site_total{i}.Coefficients{4,4};
end

% FDR correction
[h crit_p adj_ci adj_pvals_ari_total] = fdr_bh(uncorrected_pvals);

% Display significant
ROI_list(h)

% ---- Using ROIs COMBINED ------

% MAPS DATA
% Fit general linear model
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(maps_thickness{:,rois},2);
    age = maps_info{:,'age'};
    sex = maps_info{:,'sex'};
    ari = maps_info{:,'ari'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,ari);
    lm_age_sex_ari_maps_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+ari']);
    irr_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{4,4};
    irr_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{4,1};
    age_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{2,1};
    sex_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{3,1};
    age_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{2,4};
    sex_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{3,4};
    R2_total_model(f,1) = lm_age_sex_ari_maps_bigROI{f}.Rsquared.Ordinary;
end

% FDR correction
[irr_significant crit_p adj_ci irr_FDRadj_pvals] = fdr_bh(irr_uncorrected_pvals);
[age_significant crit_p adj_ci age_FDRadj_pvals] = fdr_bh(age_uncorrected_pvals);
[sex_significant crit_p adj_ci sex_FDRadj_pvals] = fdr_bh(sex_uncorrected_pvals);

results_maps_ari_bigROI = table(fields,R2_total_model,irr_significant,irr_coeffs,irr_uncorrected_pvals,irr_FDRadj_pvals,sex_significant,sex_coeffs,sex_uncorrected_pvals,sex_FDRadj_pvals,age_significant,age_coeffs,age_uncorrected_pvals,age_FDRadj_pvals);

% Display significant
fields(h)

% TEMPER LOSS
% Fit general linearl model
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(maps_thickness{:,rois},2);
    age = maps_info{:,'age'};
    sex = maps_info{:,'sex'};
    temper_loss = maps_info{:,'temper_loss'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,temper_loss);
    lm_age_sex_temper_maps_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+temper_loss']);
    uncorrected_pvals(f) = lm_age_sex_temper_maps_bigROI{f}.Coefficients{4,4};
end

% FDR correction
[h crit_p adj_ci adj_pvals_temper_maps] = fdr_bh(uncorrected_pvals);

% Display significant
fields(h)

% BDS DATA

% Fit general linear model
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(bds_thickness{:,rois},2);
    age = bds_info{:,'MRIage'};
    sex = bds_info{:,'sex'};
    ari = bds_info{:,'ARIcat'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,ari);
    lm_age_sex_ari_bds_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+ari']);
    irr_uncorrected_pvals(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{4,4};
    irr_coeffs(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{4,1};
    age_coeffs(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{2,1};
    sex_coeffs(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{3,1};
    age_uncorrected_pvals(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{2,4};
    sex_uncorrected_pvals(f,1) = lm_age_sex_ari_bds_bigROI{f}.Coefficients{3,4};
    R2_total_model(f,1) = lm_age_sex_ari_bds_bigROI{f}.Rsquared.Ordinary;
end

% FDR correction
[irr_significant crit_p adj_ci irr_FDRadj_pvals] = fdr_bh(irr_uncorrected_pvals);
[age_significant crit_p adj_ci age_FDRadj_pvals] = fdr_bh(age_uncorrected_pvals);
[sex_significant crit_p adj_ci sex_FDRadj_pvals] = fdr_bh(sex_uncorrected_pvals);

results_bds_ari_bigROI = table(fields,R2_total_model,irr_significant,irr_coeffs,irr_uncorrected_pvals,irr_FDRadj_pvals,sex_significant,sex_coeffs,sex_uncorrected_pvals,sex_FDRadj_pvals,age_significant,age_coeffs,age_uncorrected_pvals,age_FDRadj_pvals);

% Display significant
fields(h) 

% ----- Correct for CortexVol

% MAPS DATA
% Fit general linear model
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(maps_thickness{:,rois},2);
    age = maps_info{:,'age'};
    sex = maps_info{:,'sex'};
    ari = maps_info{:,'ari'};
    vol = maps_info{:,'CortexVol'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,ari,vol);
    lm_age_sex_ari_maps_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+ari+vol']);
    irr_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{4,4};
    irr_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{4,1};
    age_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{2,1};
    sex_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{3,1};
    vol_coeffs(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{5,1};
    age_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{2,4};
    sex_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{3,4};
    vol_uncorrected_pvals(f,1) = lm_age_sex_ari_maps_bigROI{f}.Coefficients{5,4};
    R2_total_model(f,1) = lm_age_sex_ari_maps_bigROI{f}.Rsquared.Ordinary;
end

% FDR correction
[irr_significant crit_p adj_ci irr_FDRadj_pvals] = fdr_bh(irr_uncorrected_pvals);
[age_significant crit_p adj_ci age_FDRadj_pvals] = fdr_bh(age_uncorrected_pvals);
[sex_significant crit_p adj_ci sex_FDRadj_pvals] = fdr_bh(sex_uncorrected_pvals);
[vol_significant crit_p adj_ci vol_FDRadj_pvals] = fdr_bh(vol_uncorrected_pvals);
results_maps_ari_bigROI = table(fields,R2_total_model,irr_significant,irr_coeffs,irr_uncorrected_pvals,irr_FDRadj_pvals,sex_significant,sex_coeffs,sex_uncorrected_pvals,sex_FDRadj_pvals,age_significant,age_coeffs,age_uncorrected_pvals,age_FDRadj_pvals,vol_significant,vol_coeffs,vol_uncorrected_pvals,vol_FDRadj_pvals);

% TEMPER LOSS
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(maps_thickness{:,rois},2);
    age = maps_info{:,'age'};
    sex = maps_info{:,'sex'};
    temper_loss = maps_info{:,'temper_loss'};
    vol = maps_info{:,'CortexVol'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,temper_loss,vol);
    lm_age_sex_temper_vol_maps_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+temper_loss+vol']);
    irr_uncorrected_pvals(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{4,4};
    irr_coeffs(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{4,1};
    age_coeffs(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{2,1};
    sex_coeffs(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{3,1};
    vol_coeffs(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{5,1};
    age_uncorrected_pvals(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{2,4};
    sex_uncorrected_pvals(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{3,4};
    vol_uncorrected_pvals(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Coefficients{5,4};
    R2_total_model(f,1) = lm_age_sex_temper_vol_maps_bigROI{f}.Rsquared.Ordinary;
end

% FDR correction
[irr_significant crit_p adj_ci irr_FDRadj_pvals] = fdr_bh(irr_uncorrected_pvals);
[age_significant crit_p adj_ci age_FDRadj_pvals] = fdr_bh(age_uncorrected_pvals);
[sex_significant crit_p adj_ci sex_FDRadj_pvals] = fdr_bh(sex_uncorrected_pvals);
[vol_significant crit_p adj_ci vol_FDRadj_pvals] = fdr_bh(vol_uncorrected_pvals);
results_maps_temper_bigROI = table(fields,R2_total_model,irr_significant,irr_coeffs,irr_uncorrected_pvals,irr_FDRadj_pvals,sex_significant,sex_coeffs,sex_uncorrected_pvals,sex_FDRadj_pvals,age_significant,age_coeffs,age_uncorrected_pvals,age_FDRadj_pvals,vol_significant,vol_coeffs,vol_uncorrected_pvals,vol_FDRadj_pvals);



% BDS DATA
% Fit general linear model
fields = fieldnames(ROI_groups);
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(bds_thickness{:,rois},2);
    age = bds_info{:,'MRIage'};
    sex = bds_info{:,'sex'};
    ari = bds_info{:,'ARIscale'};
    vol = bds_info{:,'CortexVol'};
    % Make a temporary table
    temp_T = table(avg_thickness,age,sex,ari,vol);
    lm_age_sex_ari_vol_bds_bigROI{f} = fitlm(temp_T,['avg_thickness~age+sex+ari+vol']);
    irr_uncorrected_pvals(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{4,4};
    irr_coeffs(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{4,1};
    age_coeffs(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{2,1};
    sex_coeffs(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{3,1};
    vol_coeffs(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{5,1};
    age_uncorrected_pvals(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{2,4};
    sex_uncorrected_pvals(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{3,4};
    vol_uncorrected_pvals(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Coefficients{5,4};
    R2_total_model(f,1) = lm_age_sex_ari_vol_bds_bigROI{f}.Rsquared.Ordinary;
end

% FDR correction
[irr_significant crit_p adj_ci irr_FDRadj_pvals] = fdr_bh(irr_uncorrected_pvals);
[age_significant crit_p adj_ci age_FDRadj_pvals] = fdr_bh(age_uncorrected_pvals);
[sex_significant crit_p adj_ci sex_FDRadj_pvals] = fdr_bh(sex_uncorrected_pvals);
[vol_significant crit_p adj_ci vol_FDRadj_pvals] = fdr_bh(vol_uncorrected_pvals);
results_bds_ari_bigROI = table(fields,R2_total_model,irr_significant,irr_coeffs,irr_uncorrected_pvals,irr_FDRadj_pvals,sex_significant,sex_coeffs,sex_uncorrected_pvals,sex_FDRadj_pvals,age_significant,age_coeffs,age_uncorrected_pvals,age_FDRadj_pvals,vol_significant,vol_coeffs,vol_uncorrected_pvals,vol_FDRadj_pvals);

% PLOT THICKNESS by irritable
figure;
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    irr_bds = bds_info.ARIcat == categorical({'irritable'});
    irr_maps = maps_info.ari == categorical({'irritable'});
    
    avg_thickness_bds_irr = mean(mean(bds_thickness{irr_bds,rois},2));
    avg_thickness_maps_irr = mean(mean(maps_thickness{irr_maps,rois},2));
    
    avg_thickness_bds_ctrl = mean(mean(bds_thickness{~irr_bds,rois},2));
    avg_thickness_maps_ctrl = mean(mean(maps_thickness{~irr_maps,rois},2));
    
    std_thickness_bds_irr = std(mean(bds_thickness{irr_bds,rois},2));
    std_thickness_bds_ctrl = std(mean(bds_thickness{~irr_bds,rois},2));
    
    std_thickness_maps_irr = std(mean(maps_thickness{irr_maps,rois},2));
    std_thickness_maps_ctrl = std(mean(maps_thickness{~irr_maps,rois},2));
    
    subplot(2,5,f)
    
    bar([avg_thickness_maps_irr avg_thickness_maps_ctrl; avg_thickness_bds_irr avg_thickness_bds_ctrl])
    
    hold on
    errorbar([0.85 1.15 1.85 2.15], [avg_thickness_maps_irr avg_thickness_maps_ctrl avg_thickness_bds_irr avg_thickness_bds_ctrl],[std_thickness_maps_irr std_thickness_maps_ctrl std_thickness_bds_irr std_thickness_bds_ctrl],'k.')
    xticklabels({'NU','WUSTL'})
    axis([.5 2.5 min(min([avg_thickness_maps_irr avg_thickness_bds_irr; avg_thickness_maps_ctrl avg_thickness_bds_ctrl]))-.25 max(max([avg_thickness_maps_irr avg_thickness_bds_irr; avg_thickness_maps_ctrl avg_thickness_bds_ctrl]))+.25])
    ylabel('thickness')
    legend('high ARI','low ARI')
    title(fields{f})
    
figure;
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    irr_bds = bds_info.ARIcat == categorical({'irritable'});
    irr_maps = maps_info.ari == categorical({'irritable'});
    
    avg_thickness_bds_irr = mean(mean(bds_thickness{irr_bds,rois},2));
    avg_thickness_maps_irr = mean(mean(maps_thickness{irr_maps,rois},2));
    
    avg_thickness_bds_ctrl = mean(mean(bds_thickness{~irr_bds,rois},2));
    avg_thickness_maps_ctrl = mean(mean(maps_thickness{~irr_maps,rois},2));
    
    std_thickness_bds_irr = std(mean(bds_thickness{irr_bds,rois},2));
    std_thickness_bds_ctrl = std(mean(bds_thickness{~irr_bds,rois},2));
    
    std_thickness_maps_irr = std(mean(maps_thickness{irr_maps,rois},2));
    std_thickness_maps_ctrl = std(mean(maps_thickness{~irr_maps,rois},2));
    
    subplot(2,5,f)
    
    bar([avg_thickness_maps_irr avg_thickness_maps_ctrl; avg_thickness_bds_irr avg_thickness_bds_ctrl])
    
    hold on
    errorbar([0.85 1.15 1.85 2.15], [avg_thickness_maps_irr avg_thickness_maps_ctrl avg_thickness_bds_irr avg_thickness_bds_ctrl],[std_thickness_maps_irr std_thickness_maps_ctrl std_thickness_bds_irr std_thickness_bds_ctrl],'k.')
    xticklabels({'NU','WUSTL'})
    axis([.5 2.5 min(min([avg_thickness_maps_irr avg_thickness_bds_irr; avg_thickness_maps_ctrl avg_thickness_bds_ctrl]))-.25 max(max([avg_thickness_maps_irr avg_thickness_bds_irr; avg_thickness_maps_ctrl avg_thickness_bds_ctrl]))+.25])
    ylabel('thickness')
    legend('high ARI','low ARI')
    title(fields{f})
    
end


% ROI temper_loss
figure;
for f = 1:length(fields)
    %Average ROIs
    rois = getfield(ROI_groups,fields{f});
    avg_thickness = mean(maps_thickness{:,rois},2);
   
    subplot(2,5,f)
    
    plot(maps_info.temper_loss,avg_thickness,'ko')
    
    xlabel(['MAP-DB: r = ',num2str(round(corr(avg_thickness,maps_info.temper_loss),2))])
    ylabel('thickness')
        
    title(fields{f})
    
end
