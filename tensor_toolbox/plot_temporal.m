function plot_temporal(met_specAUC_before,met_specAUC_after)
% v0: plot time series 2DEPSI
% v2: use montage
    name_mets = {'Signal1','Signal2','Signal3'};
    for i_mets = 1:size(met_specAUC_before,2)
        met_specAUC{i_mets} = cat(3,met_specAUC_before{i_mets},met_specAUC_after{i_mets});
        met_cmax(i_mets) = max(met_specAUC{i_mets}(:));
    end
%     met_cmax = [1.8e9 2.1e8];
    figure,
    for i_mets = 1:size(met_specAUC_before,2)
        subplot(3,1,i_mets)
        temp1 = permute(met_specAUC{i_mets},[2,1,4,3]);
        montage(temp1, 'Size', [2 size(met_specAUC_before{i_mets},3)], 'DisplayRange', [0 met_cmax(i_mets)]);
        %montage(temp1, 'Size', [2 size(met_specAUC_before{i_mets},3)], 'DisplayRange', [0 met_cmax(3)]/2);
        axis image off, colormap(jet(256));
        title(name_mets{i_mets});
    end

%     set(subplot(2,1,1), 'Position', [0.05, 0.69, 0.92, 0.27])
%     set(subplot(2,1,2), 'Position', [0.05, 0.37, 0.92, 0.27])
    set(subplot(3,1,1), 'Position', [0.05, 0.70, 0.92, 0.27])
    set(subplot(3,1,2), 'Position', [0.05, 0.40, 0.92, 0.27])
    set(subplot(3,1,3), 'Position', [0.05, 0.10, 0.92, 0.27])


%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 0.8]);

%     figure,
%     subplot(2,2,1)
%     imagesc(met_tempAUC_before{1}');
%     colormap gray; caxis([0 met_cmax(1)]);
%     set(gca,'xtick',[],'ytick',[]);title('before');ylabel('Pyr');
%     subplot(2,2,2)
%     imagesc(met_tempAUC_after{1}');
%     colormap gray; caxis([0 met_cmax(1)]);
%     axis off;title('after');
%     subplot(2,2,3)
%     imagesc(met_tempAUC_before{2}');
%     colormap gray; caxis([0 met_cmax(2)]);
%     set(gca,'xtick',[],'ytick',[]);ylabel('Lac');
%     subplot(2,2,4)
%     imagesc(met_tempAUC_after{2}');
%     colormap gray; caxis([0 met_cmax(2)]);
%     axis off;

end