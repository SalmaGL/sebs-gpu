function FigHandle=Graph(Data,lat,lon,GraphTitle)
Data(Data==-2)=NaN;
FigHandle=surface(lon,lat,Data,'EdgeColor','none');
shading flat;
xlim([min(min(lon)) max(max(lon))]);
ylim([min(min(lat)) max(max(lat))]);
xlabel('Longitude');
ylabel('Latitude');
title(GraphTitle);
colorbar

end
