%C = dataset that is clustered
%z = cluster to plot
%clusters = vectors with cluster assignments
%mycmap = colormap for gradinet line plots

function PlotCluster(C, z, clusters, mycmap)

cmaprange = -1:1/32:1;
cmaprange = cmaprange * nanmax(nanmax(abs(C),[],1),[],2);
[n,p] = size(C(clusters == z,:));
if n~=p, clust9 = plot(0:p-1,C(clusters == z,:),'LineWidth',2);
else
    clust9 = plot(0:p-1,C(clusters == z,:)','LineWidth',2);
end
minY = cmaprange(1,1);
maxY = cmaprange(1,end);
axis([0, size(clust9(1).YData,2)-1, minY, maxY]);
drawnow;
for i = 1:length(clust9)
    a = ones(4,3); 
    for ii = 1:size(clust9(i).YData,2)
        for iii = 1:size(cmaprange,2)
            if clust9(i).YData(ii) <= cmaprange(1,iii)
                a(:,ii) = [mycmap(min(64,iii),:)' * 255;255];
                break;
            end
        end
    end
set(clust9(i).Edge, 'ColorBinding','interpolated', 'ColorData',uint8(a),'ColorType','truecoloralpha');
end
end