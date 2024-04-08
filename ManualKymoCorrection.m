function [Tracks,PlottedTracks,TrackLabels]=ManualKymoCorrection(Tracks,Image,PlottedTracks,TrackLabels)

OK=1;
while OK==1
spaceList = {'Split a track','Delete a track','Merge a track','Add a Track','Adjust kymo scaling','Turn off the tracks','Turn on the tracks','Done'}; 
[idx, ~] = listdlg('ListString', spaceList,...
'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');
% Handle response
switch idx
    case 1
          h = msgbox({'Click near a track to select it'});
          [x,y]=ginputColour(1);
         try
             close(h)
         catch
         end
         [~,dist]=cellfun(@(z) dsearchn(z,[y(1) x(1)]),Tracks);
         [~,TrackNumber] = min(dist);
         bubbulah=plot(Tracks{TrackNumber}(:,2),Tracks{TrackNumber}(:,1),'LineWidth',4);
          h = msgbox({'Click on the point where you want to split it'});
         [x,y]=ginputColour(1);
         ClosestPoint=dsearchn(Tracks{TrackNumber},[y(1) x(1)]);
         
         Tracks{end+1}=Tracks{TrackNumber}(1:ClosestPoint,:);
         Tracks{end+1}=Tracks{TrackNumber}(ClosestPoint+1:end,:);
         Tracks{TrackNumber}=[];
         Tracks=Tracks(~cellfun('isempty',Tracks));
         try
             close(h)
         catch
         end
         delete(PlottedTracks)
         delete(TrackLabels)
         delete(bubbulah)
         [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image);
    case 2
         h = msgbox({'Click near a track to select it'});
        [x,y]=ginputColour(1);
         try
             close(h)
         catch
         end
        [~,dist]=cellfun(@(z) dsearchn(z,[y(1) x(1)]),Tracks);
        [~,TrackNumber] = min(dist);
        
         blurb1=plot(Tracks{TrackNumber}(:,2),Tracks{TrackNumber}(:,1),'LineWidth',4);
        answer = questdlg('Do you want to delete this track?','Question','Yes');
        switch answer
            case 'Yes'
                Tracks{TrackNumber}=[];
           Tracks=Tracks(~cellfun('isempty',Tracks));
        
         delete(blurb1)
         delete(PlottedTracks)
         delete(TrackLabels)
         [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image);
         
            case 'No'
           delete(blurb1)                
        end
    case 3
        h = msgbox({'Click near a two tracks whose ends are to be joined, the upper one first!'});
        [x,y]=ginputColour(2);
        [~,dist]=cellfun(@(z) dsearchn(z,[y(1) x(1)]),Tracks);
        [~,TrackNumber1] = min(dist);
        [~,dist]=cellfun(@(z) dsearchn(z,[y(2) x(2)]),Tracks);
        [~,TrackNumber2] = min(dist);
        
       blurb1=plot(Tracks{TrackNumber1}(:,2),Tracks{TrackNumber1}(:,1),'LineWidth',4);
       blurb2=plot(Tracks{TrackNumber2}(:,2),Tracks{TrackNumber2}(:,1),'LineWidth',4);
  
        answer = questdlg('Do you want to merge these tracks?','Question','Yes');
        switch answer
            case 'Yes'
              try  
            InBetween = [Tracks{TrackNumber1}(end,1):Tracks{TrackNumber2}(1,1);round(interp1([Tracks{TrackNumber1}(end,1) Tracks{TrackNumber2}(1,1)],[Tracks{TrackNumber1}(end,2) Tracks{TrackNumber2}(1,2)],Tracks{TrackNumber1}(end,1):Tracks{TrackNumber2}(1,1),'linear'))]';
            NewTrack=[Tracks{TrackNumber1}(1:end-1,:);InBetween;Tracks{TrackNumber2}(2:end,:)];
            Tracks{end+1}=NewTrack;
            Tracks{TrackNumber1}=[];
            Tracks{TrackNumber2}=[];
            Tracks=Tracks(~cellfun('isempty',Tracks));
              catch
                 h = msgbox({'It seems like the end of the first selected track does not end before the start of the second selected track.','Try again!'});
              end
            delete(blurb1)
            delete(blurb2)
            delete(PlottedTracks)
         delete(TrackLabels)
         try
             close(h)
         catch
         end
         [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image);
            case 'No'
              delete(blurb1)
            delete(blurb2)  
        end
    case 4
        h = msgbox({'Click repeatedly along a putative track, start at the earliest frame and always go down for the next point, press Enter when done'});
        

        
        Track=[];
        [x,y]=ginputColour(1);
        Track(1,1)=round(y);
        Track(1,2)=round(x);
        Point=2;
        Boink=1;
       try
        while Boink==1
            
            [x,y]=ginputColour(1);
            isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
            if isKeyPressed
               break
            end
            
            Track(Point,1)=round(y);
            Track(Point,2)=round(x);
            if exist('blorb','var')
                delete(blorb)
            end
            blorb=plot(Track(:,2),Track(:,1),'g');
            Point=Point+1;
        end
       catch
       end
       
       try
        TrackInterpolated=[Track(1,1):Track(end,1);round(interp1(Track(:,1),Track(:,2),Track(1,1):Track(end,1)))]';
        Tracks{end+1}=TrackInterpolated;
       catch 
           h=msgbox({'It seems like you didn''t always go down for the next point. Try again!'});
           pause(5)
           try
               close (h)
           catch
           end
           delete(blorb)
       end
        delete(PlottedTracks)
         delete(TrackLabels)
         delete(blorb)
         [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image);
       
           
     case 5

         Limits=get(gca,'CLim');
         prompt = {'Lower limit:','Upper limit:'};
         dlgtitle = 'Input';
         dims = [1 35];
         definput = {num2str(Limits(1)),num2str(Limits(2))};
         answer = inputdlg(prompt,dlgtitle,dims,definput);
         set(gca,'CLim',[str2num(answer{1}) str2num(answer{2})])
         
    case 6
        delete(PlottedTracks)
         delete(TrackLabels)
    case 7
        [PlottedTracks,TrackLabels]=PlotAllTracks(Tracks,Image);
            
    case 8
        OK=0;
end
end
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf,'-dpdf','Manually corrected Tracks')

save('StoredData.mat','Tracks','-append')
end