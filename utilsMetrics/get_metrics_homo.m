function metrics = get_metrics_homo(img, mask_homo, method)
% function metrics = get_metrics_homo(img, mask_homo, gt_homo, method)

    reg_homo = img(mask_homo);

    metrics.method        = method;

    % Mean and std
    metrics.mean_homo    = mean(reg_homo, 'omitnan');
    metrics.std_homo     = std(reg_homo, 'omitnan');
    metrics.cv_homo      = 100*abs(metrics.std_homo / metrics.mean_homo);
    
end