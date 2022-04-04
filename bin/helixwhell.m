function [F1]=helixwhell(w,predictor,group,savep,setting,folder_path)
if nargin <1
    
predictor=used_predictors{1};

savep=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data');
end
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
 FILENAME_whell=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_',char(setting.func{setting.helixwheel_selected_fun}),'_',char(predictor),'_whell.tif');
%%settings for testing 


%%color gradient and aa predictions 
%setting.helix_number=2
font_size=11;
for runnom=1:3
 close all
 figname=group{runnom};
 %figname=strcat(predictor,"  ",group{run});
    nametosave=strcat(savep,filesep,typeof,"_",char(setting.func{setting.helixwheel_selected_fun}),'_',predictor,'_',group{runnom},'_helixwhell.tif');
   %lines color transformation
   c1 = [0.85, 0.85, 0.85];
   c2 = [0.7, 0.7, 0.7];%/255;
   colors_p = [linspace(c1(1),c2(1),8)', linspace(c1(2),c2(2),8)', linspace(c1(3),c2(3),8)'];
   % propertis predection
   prop=[setting.func{setting.helixwheel_selected_fun}(setting.aa)]';
   if setting.helixwheel_sorting==1
      [B12,i11]=sort(prop);
      wi{1,1}= w{runnom,1}(i11);
      wi{1,2}= w{runnom,2}(i11,:);
      [B11,~]=unique(B12);
      
     
      
      %{
      valvec=B11;
      labvec=aa(i11(i12));
      ranges=max(B11)-min(B11);
      [~,pik]=max(diff(B11))
      
      for lpp=1:length(B11)
          if lpp<=pik
                 newprop(lpp,:)=colors_p1(lpp,:);
          elseif lpp>pik
              rn=length(B11)-pik
                 newprop(lpp,:)=colors_p1(lpp,:);

          end
      end
      %}
      c11 = [1, 0, 0];
      c21 = [0, 1, 0];%/255;
      colors_p1 = [linspace(c11(1),c21(1),length(B11))', linspace(c11(2),c21(2),length(B11))', linspace(c11(3),c21(3),length(B11))'];
      for qq=i11
      newprop(qq,:)=colors_p1(B11==B12(qq),:);
      end
   elseif setting.helixwheel_sorting==0
      %[B12,i11]=sort(prop);
      wi{1,1}= w{runnom,1};
      wi{1,2}= w{runnom,2};
      [B11,~]=unique(sort(prop));
      c11 = [1, 0, 0];
      c21 = [0, 1, 0];%/255;
      colors_p1 = [linspace(c11(1),c21(1),length(B11))', linspace(c11(2),c21(2),length(B11))', linspace(c11(3),c21(3),length(B11))'];
      for qq=1:20
      newprop(qq,:)=colors_p1(B11==prop(qq),:);
      end
   end
      CS=newprop;
%%line ploting for helix level 1
f1=figure('Position',[20,20,800,800],'Visible','off');
if runnom==1&&setting.helix_number>1
   title(strrep(predictor,'_','-'),'FontSize',25)
elseif runnom==1&&setting.helix_number==1
    title({strrep(predictor,'_','-');figname},'FontSize',25)
elseif runnom==2
    title(typeof,'FontSize',25)
elseif runnom==3
    title(strrep(char(setting.func{setting.helixwheel_selected_fun}),'_','-'),'FontSize',25)
%{    
elseif run==3
    axd=gca();
    axd.Colormap = CS;
    cb=colorbar('LOCATION','South');
    caxis([min(valvec) max(valvec)])
    set(cb,'XTick',valvec)
    set(cb,'XTickLabel',labvec)

%}
end
minp=0.1;
xlim([-1 1])
ylim([-1 1])
axis off
 for st=1:8
    if st<4
        [x1,y1] = pol2cart(deg2rad((0+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((0+(90*(st+1)))),0.60);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st==4
        [x1,y1] = pol2cart(deg2rad((0+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((45+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st>=5&&st<8
        [x1,y1] = pol2cart(deg2rad((45+(90*st))),0.70);
        [x2,y2] = pol2cart(deg2rad((45+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    end
end
%%pie ploting helix level 1
for st=1:4
    aanum=st;
    [x,y] = pol2cart(deg2rad((0+(90*st))),0.50);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',7)
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st)>minp,st),wi{1,1}(wi{1,2}(:,st)>minp));
    ax.Colormap = CS(wi{1,2}(:,st)>minp,:);
    set(findobj(h1,'type','text'),'fontsize',font_size);
    text(0, 1.7, string(st),'FontSize',14)
    box off
end
for st=5:8
    aanum=st;
    [x,y] = pol2cart(deg2rad((45+(90*st))),0.55);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',6)
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st)>minp,st),wi{1,1}(wi{1,2}(:,st)>minp));
    set(findobj(h1,'type','text'),'fontsize',font_size);
    ax.Colormap = CS(wi{1,2}(:,st)>minp,:);
    text(0, 1.7, string(st),'FontSize',14)
    box off
end
%% line ploting for helix level 2
if setting.helix_number >1
f2=figure('Position',[20,20,800,800],'Visible','Off');
%title for the sub fig 
title(figname,'FontSize',25)
xlim([-1 1])
ylim([-1 1])
axis off
 for st=1:8
    if st<4
        [x1,y1] = pol2cart(deg2rad((45+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((45+(90*(st+1)))),0.60);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st==4
        [x1,y1] = pol2cart(deg2rad((45+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((0+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st>=5&&st<8
        [x1,y1] = pol2cart(deg2rad((0+(90*st))),0.70);
        [x2,y2] = pol2cart(deg2rad((0+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    end
end
%%pie ploting helix level 2

for st=1:4
    aanum=st+8;
    [x,y] = pol2cart(deg2rad((45+(90*st))),0.50);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',6)
    
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st+8)>minp,st+8),wi{1,1}(wi{1,2}(:,st+8)>minp));
    set(findobj(h1,'type','text'),'fontsize',font_size);
    ax.Colormap = CS(wi{1,2}(:,st+8)>minp,:);
    text(0, 1.7, string(st+8),'FontSize',14)
    box off
end
for st=5:8
    aanum=st+8;
    [x,y] = pol2cart(deg2rad((0+(90*st))),0.55);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',6)
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st+8)>minp,st+8),wi{1,1}(wi{1,2}(:,st+8)>minp));
    set(findobj(h1,'type','text'),'fontsize',font_size);
    ax.Colormap = CS(wi{1,2}(:,st+8)>minp,:);
    text(0, 1.7, string(st+8),'FontSize',14)
    box off
end
end
%% line ploting for helix level 3
if setting.helix_number >2
f3=figure('Position',[20,20,800,800],'Visible','Off');
xlim([-1 1])
ylim([-1 1])
axis off
 for st=1:8
    if st<4
        [x1,y1] = pol2cart(deg2rad((0+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((0+(90*(st+1)))),0.60);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st==4
        [x1,y1] = pol2cart(deg2rad((0+(90*st))),0.60);
        [x2,y2] = pol2cart(deg2rad((45+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    elseif st>=5&&st<8
        [x1,y1] = pol2cart(deg2rad((45+(90*st))),0.70);
        [x2,y2] = pol2cart(deg2rad((45+(90*(st+1)))),0.70);
        line([x1 x2],[y1 y2],'LineWidth',(st/3)+1,'Color',colors_p(st,:),'LineStyle' ,':')
    end
end
%%pie ploting helix level 3

for st=1:4
    aanum=st+16;
    [x,y] = pol2cart(deg2rad((0+(90*st))),0.50);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',6)
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st+16)>minp,st+16),wi{1,1}(wi{1,2}(:,st+16)>minp));
    set(findobj(h1,'type','text'),'fontsize',font_size);
    ax.Colormap = CS(wi{1,2}(:,st+16)>minp,:);
    text(0, 1.7, string(st+16),'FontSize',14)
    box off
end
for st=5:8
    aanum=st+16;
    [x,y] = pol2cart(deg2rad((45+(90*st))),0.55);
    axes('Position',[(x+0.9)/2 (y+0.9)/2 0.15 0.15])
    set(gcf,'DefaultTextFontSize',6)
    box on
    ax = gca(); 
    h1=pie(ax,wi{1,2}(wi{1,2}(:,st+16)>minp,st+16),wi{1,1}(wi{1,2}(:,st+16)>minp));
    set(findobj(h1,'type','text'),'fontsize',font_size);
    ax.Colormap = CS(wi{1,2}(:,st+16)>minp,:);
    text(0, 1.7, string(st+16),'FontSize',14)
    box off
end
end
%% cat images, save and close figures
if setting.helix_number ==1
F1{runnom}=imcrop(frame2im(getframe(f1)),[116.5 0.5 1300 1300]);
close(f1)
elseif setting.helix_number ==2
F1{runnom} = cat(2,imcrop(frame2im(getframe(f1)),[116.5 0.5 1300 1300]),imcrop(frame2im(getframe(f2)),[116.5 0.5 1300 1300]));
close(f1)
close(f2)
elseif setting.helix_number ==3
F1{runnom} = cat(2,imcrop(frame2im(getframe(f1)),[116.5 0.5 1300 1300]),imcrop(frame2im(getframe(f2)),[116.5 0.5 1300 1300]),imcrop(frame2im(getframe(f3)),[116.5 0.5 1300 1300]));
close(f1)
close(f2)
close(f3)
end
if setting.save_intermidiet_whells==1
imwrite(F1{runnom},nametosave)
end
end
finalfig=cat(1,F1{1},F1{2},F1{3});
if setting.show_figs==1
imshow(finalfig)
end
imwrite(finalfig,FILENAME_whell)
end