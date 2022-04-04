 function plotall(processed_data,nonprocessed_table_out,plot_sum,plot_signal_disterbution,plot_corolation,plot_ratio,save_figs,folder_path)
warning('off', 'MATLAB:MKDIR:DirectoryExists');
if plot_sum==1
%plot all
analyzed_out=figure('name','analyzed out');
subplot(5,2,[1,2,3,4])
list=unique([processed_data(2:end,14)])';
las=categorical(list);
res=zeros(1,length(list));
for catnum=1:length(list)
res(catnum)=sum(string(processed_data(:,14))==list{catnum});
end
bar(las,res)
title('hit disterbution')

ppl=processed_data(2:end,1);
ppl2=nonprocessed_table_out(2:end,1);

subplot(5,2,5)
int=5;
boxplot([processed_data{2:end,int}],ppl,'PlotStyle', 'traditional')
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(processed_data{1,int},'_','-'))

subplot(5,2,7)
int=6;
boxplot([processed_data{2:end,int}],ppl)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(processed_data{1,int},'_','-'))


subplot(5,2,6)
int=9;
boxplot([processed_data{2:end,int}],ppl)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(processed_data{1,int},'_','-'))

subplot(5,2,8)
int=10;
boxplot([processed_data{2:end,int}],ppl)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(processed_data{1,int},'_','-'))


subplot(5,2,[9,10])
int=13;
boxplot([processed_data{2:end,int}],ppl)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(processed_data{1,int},'_','-'))

if save_figs==1
    try 
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'plots'))
    catch
    end
    savefig(analyzed_out,strcat(folder_path,filesep,'analesis_files',filesep,'plots',filesep,'analyzed_out.fig'))
end
end
%%
if plot_signal_disterbution==1
row_data=figure('Name','row data');

subplot(4,2,1)
int=5;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))


subplot(4,2,3)
int=7;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,5)
int=6;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,7)
int=8;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,2)
int=11;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,4)
int=13;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,6)
int=12;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))

subplot(4,2,8)
int=14;
boxplot([nonprocessed_table_out{2:end,int}],ppl2)
h=findobj(gca,'tag','Outliers'); delete(h) 
title(strrep(nonprocessed_table_out{1,int},'_','-'))
if save_figs==1
     try 
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'plots'))
    catch
    end
    savefig(row_data,strcat(folder_path,filesep,'analesis_files',filesep,'plots',filesep,'row_data.fig'))
end
end
%%
if plot_corolation==1
corelation=figure('Name','corelation');
subplot(2,2,1)
int1=5;
int2=6;

 x1=[nonprocessed_table_out{2:end,int1}];
 y1=[nonprocessed_table_out{2:end,int2}];

scatter(x1,y1,5,'filled')
 xlabel(nonprocessed_table_out{1,int1})
 ylabel(nonprocessed_table_out{1,int2})
 
subplot(2,2,3)
int1=7;
int2=8;

 x2=[nonprocessed_table_out{2:end,int1}];
 y2=[nonprocessed_table_out{2:end,int2}];

scatter(x2,y2,5,'filled')
 xlabel(nonprocessed_table_out{1,int1})
 ylabel(nonprocessed_table_out{1,int2})
 
 subplot(2,2,2)
int1=11;
int2=12;

 x3=[nonprocessed_table_out{2:end,int1}];
 y3=[nonprocessed_table_out{2:end,int2}];

scatter(x3,y3,5,'filled')
 xlabel(nonprocessed_table_out{1,int1})
 ylabel(nonprocessed_table_out{1,int2})
 
subplot(2,2,4)
int1=13;
int2=14;

 x4=[nonprocessed_table_out{2:end,int1}];
 y4=[nonprocessed_table_out{2:end,int2}];

scatter(x4,y4,5,'filled')
 xlabel(nonprocessed_table_out{1,int1})
 ylabel(nonprocessed_table_out{1,int2})
 if save_figs==1
    try 
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'plots'))
    catch
    end
    savefig(corelation,strcat(folder_path,filesep,'analesis_files',filesep,'plots',filesep,'corelation.fig'))
end
end
if plot_ratio==1
ratio=figure('Name','ratio');

int1=7;
int2=11;
names=[processed_data{2:end,3}];
exclod=strcmp([processed_data(2:end,14)],'below TH')';
xin=[processed_data{2:end,int1}];
yin=[processed_data{2:end,int2}];
%'percentiles', [0 99.9]
[~,xout] = rmoutliers(xin,'gesd','MaxNumOutliers',2);
[~,yout] = rmoutliers(yin,'gesd','MaxNumOutliers',2);
ind=logical(~xout.*~yout.*~exclod);
Xdata=xin(ind);
Ydata=yin(ind);
namesus=names(ind);
all=[Xdata Ydata];
lsto=std(all)*2;

datatype=Xdata./Ydata;
high=(datatype>=2);
low=(datatype<=(1/2));
same=((1/2)<=datatype&datatype<=2);

ef=1.1;
% discard any zero elements

allZeros = ((Xdata == 0) | (Ydata == 0));
allNegative = ((Xdata < 0) | (Ydata < 0));
upperBound = ef*max(max(Xdata(:)),max(Ydata(:)));
goodVals = ~(allZeros|allNegative);
lowerBound = 0;
hPlot = scatter(Xdata(high),Ydata(high),'.','r');
hold on
hPlot = scatter(Xdata(low),Ydata(low),'.','g');
hPlot = scatter(Xdata(same),Ydata(same),'.','k');
hold off
hAxis = get(hPlot,'parent');
lsto=hAxis.XTick(2);
lssttikto=hAxis.XTick(end);

set(hAxis,'Xlim',[lowerBound,upperBound],'Ylim',[lowerBound,upperBound]);

 line( [lowerBound,lssttikto],[lowerBound,lssttikto],'color','k','linestyle','-.','linewidth',3)

%make y line o2
 line( [0,(upperBound)/2],[0,upperBound],'color','k','linestyle','-.','linewidth',2)
 %make x lineo2
 line( [0,upperBound],[0,(upperBound)/2],'color','k','linestyle','-.','linewidth',2)
 n1=split(processed_data{1,int1},'_');
 n2=split(processed_data{1,int2},'_');
 legend({n1{1},n2{1},'same'});
 legend('boxoff');
  legend('Location','northwest');

 xlabel(strrep(processed_data{1,int1},'_','-'))
 ylabel(strrep(processed_data{1,int2},'_','-'))
dcm_obj = datacursormode(ratio);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,namesus,Xdata,Ydata})


 if save_figs==1
    try 
    mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'plots'))
    catch
    end
    savefig(ratio,strcat(folder_path,filesep,'analesis_files',filesep,'plots',filesep,'ratio.fig'))
end
end

 end