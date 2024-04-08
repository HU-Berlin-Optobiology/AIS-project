function [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image)


for tt=1:length(Tracks)
    eval(['PlottedTracks(' num2str(tt) ')=plot(Tracks{tt}(:,2),Tracks{tt}(:,1),''LineWidth'',2);'])
    PlotColor=get(PlottedTracks(tt),'Color');
    eval(['TrackLabels(' num2str(tt) ')=text(Tracks{tt}(round(length(Tracks{tt})/2),2),Tracks{tt}(round(length(Tracks{tt})/2),1),num2str(tt),''Color'',PlotColor);'])
end
xlabel('Distance [pixels]')
ylabel('Frame [#]')
title('Tracks')

xlim([1 size(Image.data,2)])
set(gca,'YDir',"reverse")
grid minor
xlabel('Dist')
ylabel('Frame')
%%%%