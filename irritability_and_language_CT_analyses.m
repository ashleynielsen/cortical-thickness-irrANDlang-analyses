% Data loaded from MAPS_irr_lang.txt
% load maps_data.mat
load('/projects/p30906/data/maps_bds/maps/maps_data.mat')
load('/projects/p30906/data/maps_bds/bds/bds_data_02.17.20.mat')

warning off
% Want to grab the ids for those that are high/low on language and
% irritability

%%%% MAPS %%%%
% fit a linear model
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
maps_singleROI_fullModel_temper_loss = full_regression(more_maps_info, more_maps_info(:,'temper_loss'),maps_thickness, 1);
% categorize papa_ari
temp_papa_ari = table(categorical(more_maps_info{:,'papa_ari'}>ceil(nanmedian(more_maps_info.papa_ari))),'VariableNames',{'papa_ari'});
maps_singleROI_fullModel_papa_ari = full_regression(more_maps_info, temp_papa_ari, maps_thickness, 1);
maps_singleROI_fullModel_ex_vocab = full_regression(more_maps_info, more_maps_info(:,'ex_vocab'),maps_thickness, 1);
maps_singleROI_fullModel_sent_struct = full_regression(more_maps_info, more_maps_info(:,'sent_struct'),maps_thickness, 1);
maps_singleROI_fullModel_follow_dir = full_regression(more_maps_info, more_maps_info(:,'follow_dir'),maps_thickness, 1);


% regress out some things
maps_thickness_regress = remove_regressors(more_maps_info,maps_thickness);

% temper loss from MAP-DB
[maps_hi_temper_loss_ids, maps_lo_temper_loss_ids, maps_meanCT_hi_temper_loss, maps_meanCT_lo_temper_loss, ...
    maps_meanCT_hiVSlo_temper_loss] = calculate_hi_lo_diff(more_maps_info, more_maps_info{:,'temper_loss'}, ...
    maps_thickness_regress);

% papa - ari
[maps_hi_papa_ari_ids, maps_lo_papa_ari_ids, maps_meanCT_hi_papa_ari, maps_meanCT_lo_papa_ari, ...
    maps_meanCT_hiVSlo_papa_ari] = calculate_hi_lo_diff(more_maps_info, more_maps_info{:,'papa_ari'}, ...
    maps_thickness_regress);

% expressive vocabulary
[maps_hi_ex_vocab_ids, maps_lo_ex_vocab_ids, maps_meanCT_hi_ex_vocab, maps_meanCT_lo_ex_vocab, ...
    maps_meanCT_hiVSlo_ex_vocab] = calculate_hi_lo_diff(more_maps_info, more_maps_info{:,'ex_vocab'}, ...
    maps_thickness_regress);

% sentence structure
[maps_hi_sent_struct_ids, maps_lo_sent_struct_ids, maps_meanCT_hi_sent_struct, maps_meanCT_lo_sent_struct, ...
    maps_meanCT_hiVSlo_sent_struct] = calculate_hi_lo_diff(more_maps_info, more_maps_info{:,'sent_struct'}, ...
    maps_thickness_regress);

% follow directions
[maps_hi_follow_dir_ids, maps_lo_follow_dir_ids, maps_meanCT_hi_follow_dir, maps_meanCT_lo_follow_dir, ...
    maps_meanCT_hiVSlo_follow_dir] = calculate_hi_lo_diff(more_maps_info, more_maps_info{:,'follow_dir'}, ...
    maps_thickness_regress);

%%% BDS %%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Full linear Model
temp_papa_ari = table(categorical(more_bds_info{:,'papa_ari'}>ceil(nanmedian(more_bds_info.papa_ari))),'VariableNames',{'papa_ari'});
bds_singleROI_fullModel_papa_ari = full_regression(more_bds_info, temp_papa_ari, bds_thickness, 1);
bds_singleROI_fullModel_verbal = full_regression(more_bds_info, more_bds_info(:,'kb2_vss'),bds_thickness, 1);


%
results = full_irr_lang_regression(more_bds_info,temp_papa_ari, more_bds_info(:,'kb2_vss'),bds_thickness, 1);

% Regress out some stuff
bds_thickness_regress = remove_regressors(more_bds_info,bds_thickness);
    % Var1 - 

% papa - ari
[bds_hi_papa_ari_ids, bds_lo_papa_ari_ids, bds_meanCT_hi_papa_ari, bds_meanCT_lo_papa_ari, ...
    bds_meanCT_hiVSlo_papa_ari] = calculate_hi_lo_diff( more_bds_info, more_bds_info{:,'papa_ari'}, ...
    bds_thickness_regress);

% KBIT2 - verbal standard scale
[bds_hi_verbal_ids, bds_lo_verbal_ids, bds_meanCT_hi_verbal, bds_meanCT_lo_verbal, ...
    bds_meanCT_hiVSlo_verbal] = calculate_hi_lo_diff( more_bds_info, more_bds_info{:,'kb2_vss'}, ...
    bds_thickness_regress);

% REPORT


% Correlated patterns across sites and variables
var_corr = corr([maps_meanCT_hiVSlo_temper_loss' maps_meanCT_hiVSlo_papa_ari' -maps_meanCT_hiVSlo_ex_vocab' -maps_meanCT_hiVSlo_sent_struct' -maps_meanCT_hiVSlo_follow_dir' bds_meanCT_hiVSlo_papa_ari' -bds_meanCT_hiVSlo_verbal']);
figure; imagesc(var_corr)
colorbar
caxis([-.2 .8])
xticklabels([more_maps_info.Properties.VariableNames(3:7),more_bds_info.Properties.VariableNames(5:6)])
xtickangle(45)
yticklabels([more_maps_info.Properties.VariableNames(3:7),more_bds_info.Properties.VariableNames(5:6)])


% % Find the ROIS that differ the most
% % temper_loss
% [blah, temp_idx] = sort(abs(meanCT_hiVSlo_temper_loss),'descend');
% temper_loss_rois = temp_idx(1:15);
% 
% % papa ari
% [blah, temp_idx] = sort(abs(meanCT_hiVSlo_papa_ari),'descend');
% papa_ari_rois = temp_idx(1:15);
% 
% % expressive vocab
% [blah, temp_idx] = sort(abs(meanCT_hiVSlo_ex_vocab),'descend');
% ex_vocab_rois = temp_idx(1:15);
% 
% % sentence structure
% [blah, temp_idx] = sort(abs(meanCT_hiVSlo_sent_struct),'descend');
% sent_struct_rois = temp_idx(1:15);
% 
% % follow_dir
% [blah, temp_idx] = sort(abs(meanCT_hiVSlo_follow_dir),'descend');
% follow_dir_rois = temp_idx(1:15);
% 
% % find overlap
% bestROI_overlap(1,2) = numel(intersect(temper_loss_rois,papa_ari_rois));
% bestROI_overlap(1,3) = numel(intersect(temper_loss_rois,ex_vocab_rois));
% bestROI_overlap(1,4) = numel(intersect(temper_loss_rois,sent_struct_rois));
% bestROI_overlap(1,5) = numel(intersect(temper_loss_rois,follow_dir_rois));
% bestROI_overlap(2,1) = numel(intersect(temper_loss_rois,papa_ari_rois));
% bestROI_overlap(3,1) = numel(intersect(temper_loss_rois,ex_vocab_rois));
% bestROI_overlap(4,1) = numel(intersect(temper_loss_rois,sent_struct_rois));
% bestROI_overlap(5,1) = numel(intersect(temper_loss_rois,follow_dir_rois));
% 
% bestROI_overlap(2,3) = numel(intersect(papa_ari_rois,ex_vocab_rois));
% bestROI_overlap(2,4) = numel(intersect(papa_ari_rois,sent_struct_rois));
% bestROI_overlap(2,5) = numel(intersect(papa_ari_rois,follow_dir_rois));
% bestROI_overlap(3,2) = numel(intersect(papa_ari_rois,ex_vocab_rois));
% bestROI_overlap(4,2) = numel(intersect(papa_ari_rois,sent_struct_rois));
% bestROI_overlap(5,2) = numel(intersect(papa_ari_rois,follow_dir_rois));
% 
% bestROI_overlap(3,4) = numel(intersect(ex_vocab_rois,sent_struct_rois));
% bestROI_overlap(3,5) = numel(intersect(ex_vocab_rois,follow_dir_rois));
% bestROI_overlap(4,3) = numel(intersect(ex_vocab_rois,sent_struct_rois));
% bestROI_overlap(5,3) = numel(intersect(ex_vocab_rois,follow_dir_rois));
% 
% bestROI_overlap(4,5) = numel(intersect(sent_struct_rois,follow_dir_rois));
% bestROI_overlap(5,4) = numel(intersect(sent_struct_rois,follow_dir_rois));

% First we need to find the median
function [hi_ids, lo_ids, meanCT_hi, meanCT_lo, meanCT_hiVSlo] = calculate_hi_lo_diff(T_info,T_voi,T_CT)
hi_ids = T_info.Properties.RowNames(T_voi > ceil(nanmedian(T_voi)));
lo_ids = T_info.Properties.RowNames(T_voi <= ceil(nanmedian(T_voi)));

meanCT_hi = mean(T_CT{hi_ids,:});
meanCT_lo = mean(T_CT{lo_ids,:});
meanCT_hiVSlo = meanCT_hi - meanCT_lo;
end

function residuals = remove_regressors(T_info, T_CT)
residuals = T_CT;
for n = 1:size(T_CT,2)
    y = T_CT{:,n};
    tbl = table(y, T_info.age, categorical(T_info.sex), T_info.time_btw, T_info.T1qc, T_info.cortex_vol);
    mdl = fitlm(tbl,'y~Var2+Var3+Var4+Var5+Var6');
    residuals{:,n} = mdl.Residuals.Raw;
end
end

function results = full_regression(T_info,T_voi, T_CT,report_on)
voi = T_voi{:,:};
uncorrected_pvals = ones(size(T_CT,2),1);
coefficient = ones(size(T_CT,2),1);
total_model_r2 = ones(size(T_CT,2),1);
for n = 1:size(T_CT,2)
    y = T_CT{:,n};
    tbl = table(y, voi, T_info.age, categorical(T_info.sex), T_info.time_btw, T_info.T1qc, T_info.cortex_vol );
    mdl = fitlm(tbl,'y~voi+Var3+Var4+Var5+Var6+Var7');
    uncorrected_pvals(n) = mdl.Coefficients{2,4};
    coefficient(n) = mdl.Coefficients{2,1};
    total_model_r2(n) = mdl.Rsquared.Ordinary;
end
results = table(coefficient,uncorrected_pvals,total_model_r2);
results.Properties.RowNames = T_CT.Properties.VariableNames';

% creat report
if report_on
    disp(['Report for: ',T_voi.Properties.VariableNames{1,:}])
    disp('------------------------------------------------------------------------------------')
    disp('ROIs w/ p<0.05:')
    sig_rois = results.Properties.RowNames(results.uncorrected_pvals<0.05);
    disp(['#',9,'ROI Name',9,9,9,9,9,'p-value',9,9,'Coefficient'])
    for n = 1:length(sig_rois)
        num_tabs = 9.*ones(1,6-floor(length(sig_rois{n})/8));
        disp([num2str(n),9,sig_rois{n},num_tabs,num2str(round(results{sig_rois{n},'uncorrected_pvals'},4)),9,9,num2str(round(results{sig_rois{n},'coefficient'},3))])
    end
    disp('-------------------------------------------------------------------------------------')
end
end


function results = full_irr_lang_regression(T_info,T_voi_i,T_voi_l, T_CT,report_on)
voi_l = T_voi_l{:,:};
voi_i = T_voi_i{:,:};

uncorrected_pvals = ones(size(T_CT,2),1);
coefficient = ones(size(T_CT,2),1);
total_model_r2 = ones(size(T_CT,2),1);
for n = 1:size(T_CT,2)
    y = T_CT{:,n};
    tbl = table(voi_i,voi_l, T_info.age, categorical(T_info.sex), T_info.time_btw, T_info.T1qc, T_info.cortex_vol,y,'VariableNames',{'irr','lang','age','sex','time','T1qc','vol','CT'}); 
    mdl = stepwiselm(tbl,'interactions');
    uncorrected_pvals(n) = mdl.Coefficients{2,4};
    coefficient(n) = mdl.Coefficients{2,1};
    total_model_r2(n) = mdl.Rsquared.Ordinary;
end
results = table(coefficient,uncorrected_pvals,total_model_r2);
results.Properties.RowNames = T_CT.Properties.VariableNames';

% creat report
if report_on
    disp(['Report for: ',T_voi.Properties.VariableNames{1,:}])
    disp('------------------------------------------------------------------------------------')
    disp('ROIs w/ p<0.05:')
    sig_rois = results.Properties.RowNames(results.uncorrected_pvals<0.05);
    disp(['#',9,'ROI Name',9,9,9,9,9,'p-value',9,9,'Coefficient'])
    for n = 1:length(sig_rois)
        num_tabs = 9.*ones(1,6-floor(length(sig_rois{n})/8));
        disp([num2str(n),9,sig_rois{n},num_tabs,num2str(round(results{sig_rois{n},'uncorrected_pvals'},4)),9,9,num2str(round(results{sig_rois{n},'coefficient'},3))])
    end
    disp('-------------------------------------------------------------------------------------')
end
end

