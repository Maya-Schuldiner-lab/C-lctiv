function plotinganotation(hfig1fraim,hfig2fraim,hfig3fraim,label,group1,group2,hitlist,FILENAME_anotated,setting,hitlistsize)
catim=cat(1,hfig1fraim,hfig2fraim,hfig3fraim);
if setting.show_figs==1
figure
else
figure('Visible','Off')
end
himage=imshow(catim);
xmin=[0.243734939759033 0.243734939759033] ;
xmax=0.719136546184729;
xdel=(xmax-xmin(1))/(setting.set_TMD_size-1);
%xdel=0.0190763052208833;
ydel=0.28306342780027;
%ymax=0.953461538461538;
ymin=[0.387334682860997 0.367334682860997];
                                if setting.get_TMD==1&&setting.get_SS==0&&strcmp(string(setting.TMD_to_check),'maxim')
                        typeof='lest TMD';            
                    elseif setting.get_TMD==0&&setting.get_SS==1
                        typeof='signal sequnce';
                    elseif setting.get_TMD==1&&setting.get_SS==0&&setting.TMD_to_check==1
                        typeof='first TMD';
                    elseif setting.get_TMD==1&&setting.get_SS==0&&setting.TMD_to_check==2
                        typeof='second TMD';
                        elseif setting.get_TMD==1&&setting.get_SS==0&&setting.TMD_to_check==3
                        typeof='third TMD';
                         elseif setting.get_TMD==1&&setting.get_SS==0&&setting.TMD_to_check==4
                        typeof='fourth TMD';
                        elseif setting.get_cyto==1&&setting.TMD_to_check==1
                        typeof='first cyto region';
                        elseif setting.get_cyto==1&&setting.TMD_to_check==2
                        typeof='second cyto region';
                        elseif setting.get_cyto==1&&setting.TMD_to_check==3
                        typeof='third cyto region';
                        elseif setting.get_luminal==1&&setting.TMD_to_check==1
                        typeof='first luminal region';
                        elseif setting.get_luminal==1&&setting.TMD_to_check==2
                        typeof='second luminal region';
                        elseif setting.get_luminal==1&&setting.TMD_to_check==3
                        typeof='third luminal region';
                    elseif setting.get_TMD==1&&setting.get_SS==1&&setting.TMD_to_check==1
                        typeof='first TMD or signal sequnce';
                    end
hh=title([strrep(label,'_','-'),'----',typeof]);
hh.FontSize = 9;
annotation('textbox',[0.042164658634538,0.814690026279409,0.126506026514084,0.045148248653206],'String',group1,'EdgeColor','none', 'Color','K');
annotation('textbox',[0.042164658634538,0.535714285039518,0.126506026514084,0.045148248653206],'String',group2,'EdgeColor','none', 'Color','K');
annotation('textbox',[0.042164658634538,0.248652290430353,0.126506026514084,0.045148248653206],'String','same','EdgeColor','none', 'Color','K');

if ~isempty(hitlist)
for yy=1:3
   clear length1 colors_p B i
   length1 = length(unique(hitlistsize{yy}));
   red = [1, 0, 0];
   pink = [0.9290, 0.6940, 0.1250];%/255;
   colors_p = [linspace(red(1),pink(1),length1)', linspace(red(2),pink(2),length1)', linspace(red(3),pink(3),length1)'];
   [B,~]=sort(unique(hitlistsize{yy}));
   
    for xx=1:length(hitlist{yy})
      
        annotation('textarrow',xmin+(xdel*(hitlist{yy}(xx)-1)),ymin+(ydel*(yy-1)),'TextColor','k','FontSize',5,'HorizontalAlignment','center','String',string(hitlist{yy}(xx)),'Color',colors_p(find(hitlistsize{yy}(xx)==B),:),'HeadWidth',2,'HeadLength',2);
        %annotation('arrow',xmin+(xdel*(hitlist{yy}(xx)-1)),ymin+(ydel*(yy-1)),'Color',colors_p(find(hitlistsize{yy}(xx)==B),:),'HeadWidth',2,'HeadLength',2);
    end
end
end
                saveas(himage,FILENAME_anotated)
              %  close(himage)

end