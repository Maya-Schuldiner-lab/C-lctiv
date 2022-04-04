function [processed_data,nonprocessed_table_out,group1,group2]=...
    dotblot_analisis(setting,row_data,folder_path)

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%groups names
groups=fieldnames(row_data);
group1=groups{1,1};
group2=groups{2,1};


if setting.use_madian==1
    out_name_string_p1=' median based ';
else
    out_name_string_p1=' mean based ';
end
if setting.LC_global_normalization==1
    out_name_string_p2=' global normalization ';
else
    out_name_string_p2=' plate normalization ';
end
if setting.use_abselute_deviation==1
    out_name_string_p3=' absolute deviation ';
else
    out_name_string_p3=' standard deviation ';
end

if setting.pre_normalization_type==0
    out_name_string_p4=' no normalization ';
elseif setting.pre_normalization_type==1
    out_name_string_p4=' min max normalization ';
    %elseif setting.pre_normalization_type==2
    %    out_name_string_p4=' microarray normalization ';
end

if setting.LC_normalization_factor==0
    out_name_string_p5=' 1k h3 normalization ';
else
    out_name_string_p5=' mean h3 normalization ';
end
if setting.use_dot_stat_median==0
    out_name_string_p6=' spot mean ';
    
else
    out_name_string_p6=' spot median ';
end
out_name_string=strcat(out_name_string_p1,out_name_string_p2,out_name_string_p3,out_name_string_p4,out_name_string_p5,out_name_string_p6);
%difine preset varabols, row number to ABC convertion vector, chanel
%counters and heders for the out tabels
row_number2abc_convert_vector='ABCDEFGHIGKLMNOPQRSTUVWXYZ';

plate_counter1=1;

nonprocessed_data{1,1}='plate number';
nonprocessed_data{1,2}='384 cordinations';
nonprocessed_data{1,3}=[group1 ' R1 LC'];
nonprocessed_data{1,4}=[group1 ' R1 SA'];
nonprocessed_data{1,5}=[group1 ' R2 LC'];
nonprocessed_data{1,6}=[group1 ' R2 SA'];
nonprocessed_data{1,7}=[group1 ' R3 LC'];
nonprocessed_data{1,8}=[group1 ' R3 SA'];
nonprocessed_data{1,9}=[group2 ' R1 LC'];
nonprocessed_data{1,10}=[group2 ' R1 SA'];
nonprocessed_data{1,11}=[group2 ' R2 LC'];
nonprocessed_data{1,12}=[group2 ' R2 SA'];
nonprocessed_data{1,13}=[group2 ' R3 LC'];
nonprocessed_data{1,14}=[group2 ' R3 SA'];

temp_sum_table1{1,1}='plate';
temp_sum_table1{1,2}='384 cord';
temp_sum_table1{1,3}='gene name';
temp_sum_table1{1,4}='ORF';

processed_data{1,1}='plate';
processed_data{1,2}='384 cord';
processed_data{1,3}='gene name';
processed_data{1,4}='ORF';
processed_data{1,5}=[group1 '_LC'];
processed_data{1,6}=[group1 '_SA'];
processed_data{1,7}=[group1 '_norm'];
processed_data{1,9}=[group2 '_LC'];
processed_data{1,10}=[group2 '_SA'];
processed_data{1,11}=[group2 '_norm'];
processed_data{1,13}=[group1 '/' group2];
processed_data{1,14}='hit';

%pre processing
if setting.use_dot_stat_median==1
    signal_col='median';
elseif setting.use_dot_stat_median==0
    signal_col='mean';
elseif setting.use_dot_stat_median==3
    signal_col='Median_minus_Median';
end
if setting.use_OD_read ==1
    LCreadfild='OD_Read';
else
    LCreadfild=['F680_',signal_col];
end
counter=1;
for ant=1:size(row_data.(string(groups(1))))
    if ~isempty(row_data.(string(groups(1)))(ant).plate_number)&&~isempty(row_data.(string(groups(1)))(ant).F680_median_R1)&&~isempty(row_data.(string(groups(1)))(ant).F680_median_R2)...
            &&~isempty(row_data.(string(groups(2)))(ant).F680_median_R1)&& ~isempty(row_data.(string(groups(2)))(ant).F680_median_R2)
        %% add here extra cleanup gates like snr
        %% 
        data_no_beckgrund{counter,1}=row_data.(string(groups(1)))(ant).plate_number;
        data_no_beckgrund{counter,2}=row_data.(string(groups(1)))(ant).sample_cord;
        
        %680r1 beckgrund removed
        data_no_beckgrund{counter,3}=row_data.(string(groups(1)))(ant).([LCreadfild,'_R1']);
        %800r1 beckgrund removed
        data_no_beckgrund{counter,4}=row_data.(string(groups(1)))(ant).(['F800_',signal_col,'_R1']);
        %680r2 beckgrund removed
        data_no_beckgrund{counter,5}=row_data.(string(groups(1)))(ant).([LCreadfild,'_R2']);
        %800r2 beckgrund removed
        data_no_beckgrund{counter,6}=row_data.(string(groups(1)))(ant).(['F800_',signal_col,'_R2']);
        %{
        if width(row_data.(string(groups(1))))>44
            
            %680r3 beckgrund removed
            data_no_beckgrund{counter,7}=row_data.(string(groups(1)))(ant).(['F680_',signal_col,'_R3']);
            %800r3 beckgrund removed
            data_no_beckgrund{counter,8}=row_data.(string(groups(1)))(ant).(['F800_',signal_col,'_R3']);
        end
        %}
        
        %680r1 beckgrund removed
        data_no_beckgrund{counter,9}=row_data.(string(groups(2)))(ant).([LCreadfild,'_R1']);
        %800r1 beckgrund removed
        data_no_beckgrund{counter,10}=row_data.(string(groups(2)))(ant).(['F800_',signal_col,'_R1']);
        %680r2 beckgrund removed
        data_no_beckgrund{counter,11}=row_data.(string(groups(2)))(ant).([LCreadfild,'_R2']);
        %800r2 beckgrund removed
        data_no_beckgrund{counter,12}=row_data.(string(groups(2)))(ant).(['F800_',signal_col,'_R2']);
        %{
        if width(row_data.(string(groups(1))))>44
            %680r3 beckgrund removed
            data_no_beckgrund{counter,13}=row_data.(string(groups(2)))(ant).(['F680_',signal_col,'_R3']);
            %800r3 beckgrund removed
            data_no_beckgrund{counter,14}=row_data.(string(groups(2)))(ant).(['F800_',signal_col,'_R3']);
            splist=[3,4,5,6,7,8,9,10,11,12,13,14];
        else
        %}
            data_no_beckgrund{counter,13}=[];
            data_no_beckgrund{counter,14}=[];
            splist=[3,4,5,6,9,10,11,12];
        %end
        counter=counter+1;
    end
    %textwaitbar(ant, size(row_data.(string(groups(1))),1), 'reformating')

end

if  setting.pre_normalization_type==1
    platt=unique([data_no_beckgrund{:,1}]);
    
    %prelocation
    plate_index=cell(length(platt),1);
    vforpo=cell(length(platt),max(splist));
    
    for plplp=1:length(platt)
        plate_index{plplp}=[data_no_beckgrund{:,1}]==platt(plplp);
        for sp=splist
            vforpo{plplp,sp}=[data_no_beckgrund{plate_index{plplp}',sp}];
            mins(plplp,sp)=min(vforpo{plplp,sp});
            maxs(plplp,sp)=max(vforpo{plplp,sp});
        end
    end
    for spp=splist
        minall(spp)=min(mins(:,spp));
        maxall(spp)=max(maxs(:,spp));
    end
    for wew=1:counter-1
        plo=find((data_no_beckgrund{wew,1}==platt)', 1, 'first');
        for spd=splist
            data_no_beckgrund{wew,spd}=((data_no_beckgrund{wew,spd}-mins(plo,spd))/(maxs(plo,spd)-mins(plo,spd)))*(maxall(spd)-minall(spd))+minall(spd);
        end
        textwaitbar(wew, counter-1, 'runing plate normalization')
        
    end
    
elseif setting.pre_normalization_type==0
end

disp('runing the analysis')
%prelocation
plate_list=zeros(size(data_no_beckgrund,1),1);
group1_LC_1=zeros(size(data_no_beckgrund,1),1);
group1_LC_2=zeros(size(data_no_beckgrund,1),1);
group1_SA_1=zeros(size(data_no_beckgrund,1),1);
group1_SA_2=zeros(size(data_no_beckgrund,1),1);
group2_LC_1=zeros(size(data_no_beckgrund,1),1);
group2_LC_2=zeros(size(data_no_beckgrund,1),1);
group2_SA_1=zeros(size(data_no_beckgrund,1),1);
group2_SA_2=zeros(size(data_no_beckgrund,1),1);
group1_LC_rtio=zeros(size(data_no_beckgrund,1),1);
group1_SA_rtio=zeros(size(data_no_beckgrund,1),1);
group2_LC_rtio=zeros(size(data_no_beckgrund,1),1);
group2_SA_rtio=zeros(size(data_no_beckgrund,1),1);


for read_counter=1:size(data_no_beckgrund,1)
    %readouts
    plate_list(read_counter,1)=str2double(data_no_beckgrund{read_counter,1});
    group1_LC_1(read_counter,1)=data_no_beckgrund{read_counter,3};
    group1_LC_2(read_counter,1)=data_no_beckgrund{read_counter,5};
    group1_SA_1(read_counter,1)=data_no_beckgrund{read_counter,4};
    group1_SA_2(read_counter,1)=data_no_beckgrund{read_counter,6};
    group2_LC_1(read_counter,1)=data_no_beckgrund{read_counter,9};
    group2_LC_2(read_counter,1)=data_no_beckgrund{read_counter,11};
    group2_SA_1(read_counter,1)=data_no_beckgrund{read_counter,10};
    group2_SA_2(read_counter,1)=data_no_beckgrund{read_counter,12};
    %rapets ratio
    group1_LC_rtio(read_counter,1)=data_no_beckgrund{read_counter,3}/data_no_beckgrund{read_counter,5};
    group1_SA_rtio(read_counter,1)=data_no_beckgrund{read_counter,4}/data_no_beckgrund{read_counter,6};
    group2_LC_rtio(read_counter,1)=data_no_beckgrund{read_counter,9}/data_no_beckgrund{read_counter,11};
    group2_SA_rtio(read_counter,1)=data_no_beckgrund{read_counter,10}/data_no_beckgrund{read_counter,12};
end
% calculet the mean and standart deviation for each ratio aoly for valis
% valieus
if setting.use_madian==0
    group1_LC_rtio_mean=nanmean(group1_LC_rtio(isfinite(group1_LC_rtio)));
    group1_SA_rtio_mean=nanmean(group1_SA_rtio(isfinite(group1_SA_rtio)));
    group2_LC_rtio_mean=nanmean(group2_LC_rtio(isfinite(group2_LC_rtio)));
    group2_SA_rtio_mean=nanmean(group2_SA_rtio(isfinite(group2_SA_rtio)));
elseif setting.use_madian==1
    group1_LC_rtio_mean=nanmedian(group1_LC_rtio(isfinite(group1_LC_rtio)));
    group1_SA_rtio_mean=nanmedian(group1_SA_rtio(isfinite(group1_SA_rtio)));
    group2_LC_rtio_mean=nanmedian(group2_LC_rtio(isfinite(group2_LC_rtio)));
    group2_SA_rtio_mean=nanmedian(group2_SA_rtio(isfinite(group2_SA_rtio)));
end

if setting.use_abselute_deviation==0
    group1_LC_rtio_stds=nanstd(group1_LC_rtio(isfinite(group1_LC_rtio)));
    group1_SA_rtio_stds=nanstd(group1_SA_rtio(isfinite(group1_SA_rtio)));
    group2_LC_rtio_stds=nanstd(group2_LC_rtio(isfinite(group2_LC_rtio)));
    group2_SA_rtio_stds=nanstd(group2_SA_rtio(isfinite(group2_SA_rtio)));
    
elseif setting.use_abselute_deviation==1
    group1_LC_rtio_stds=mad(group1_LC_rtio(isfinite(group1_LC_rtio)));
    group1_SA_rtio_stds=mad(group1_SA_rtio(isfinite(group1_SA_rtio)));
    group2_LC_rtio_stds=mad(group2_LC_rtio(isfinite(group2_LC_rtio)));
    group2_SA_rtio_stds=mad(group2_SA_rtio(isfinite(group2_SA_rtio)));
end

% calculate the dtandart score for each ratio and teh SA/LC ratios ratio only for valid vlaeu
%prelocate
group1_LC_STS=zeros(size(group1_LC_rtio,1),1);
group1_SA_STS=zeros(size(group1_LC_rtio,1),1);
group1_STS_rat=zeros(size(group1_LC_rtio,1),1);
group2_LC_STS=zeros(size(group1_LC_rtio,1),1);
group2_SA_STS=zeros(size(group1_LC_rtio,1),1);
group2_STS_rat=zeros(size(group1_LC_rtio,1),1);

for rep_ratio_counter=1:size(group1_LC_rtio,1)
    if isfinite(group1_LC_rtio(rep_ratio_counter))
        group1_LC_STS(rep_ratio_counter,1)=(group1_LC_rtio(rep_ratio_counter)-group1_LC_rtio_mean)/group1_LC_rtio_stds;
    end
    if isfinite(group1_SA_rtio(rep_ratio_counter))
        group1_SA_STS(rep_ratio_counter,1)=(group1_SA_rtio(rep_ratio_counter)-group1_SA_rtio_mean)/group1_SA_rtio_stds;
    end
    if isfinite(group1_LC_rtio(rep_ratio_counter))&&isfinite(group1_SA_rtio(rep_ratio_counter))
        group1_STS_rat(rep_ratio_counter,1)=group1_LC_STS(rep_ratio_counter,1)/group1_SA_STS(rep_ratio_counter,1);
    end
    if isfinite(group2_LC_rtio(rep_ratio_counter))
        group2_LC_STS(rep_ratio_counter,1)=(group2_LC_rtio(rep_ratio_counter)-group2_LC_rtio_mean)/group2_LC_rtio_stds;
    end
    if isfinite(group2_SA_rtio(rep_ratio_counter))
        group2_SA_STS(rep_ratio_counter,1)=(group2_SA_rtio(rep_ratio_counter)-group2_SA_rtio_mean)/group2_SA_rtio_stds;
    end
    if isfinite(group2_LC_rtio(rep_ratio_counter))&&isfinite(group2_SA_rtio(rep_ratio_counter))
        group2_STS_rat(rep_ratio_counter,1)=group2_LC_STS(rep_ratio_counter,1)/group2_SA_STS(rep_ratio_counter,1);
    end
    
end

if setting.use_madian==0
    group1_STS_rat_mean=nanmean(group1_STS_rat(isfinite(group1_STS_rat)));
    group2_STS_rat_mean=nanmean(group2_STS_rat(isfinite(group2_STS_rat)));
elseif  setting.use_madian==1
    group1_STS_rat_mean=nanmedian(group1_STS_rat(isfinite(group1_STS_rat)));
    group2_STS_rat_mean=nanmedian(group2_STS_rat(isfinite(group2_STS_rat)));
end

if setting.use_abselute_deviation==0
    % calculate the mean and standart deviation for the SA/LC ratios ratio
    group1_STS_rat_stds=nanstd(group1_STS_rat(isfinite(group1_STS_rat)));
    group2_STS_rat_stds=nanstd(group2_STS_rat(isfinite(group2_STS_rat)));
elseif setting.use_abselute_deviation==1
    group1_STS_rat_stds=mad(group1_STS_rat(isfinite(group1_STS_rat)));
    group2_STS_rat_stds=mad(group2_STS_rat(isfinite(group2_STS_rat)));
end


% check if between repets standart scors ratios are outliers from the
% normal disterbutions of ratios in case that all are valid numbers if both
% ratios are not sig diffrent from the normal disterbution assing 1 to the
% oklist with the correct index
%prelocation
ok_list=zeros(size(group1_LC_rtio,1),1);

for zscore_counter=1:size(group1_LC_rtio,1)
    group1_rep_zscore=(group1_STS_rat(zscore_counter)-group1_STS_rat_mean)/group1_STS_rat_stds;
    group2_rep_zscore=(group2_STS_rat(zscore_counter)-group2_STS_rat_mean)/group2_STS_rat_stds;
    if isfinite(group1_LC_rtio(rep_ratio_counter))&&isfinite(group1_SA_rtio(rep_ratio_counter))&&isfinite(group2_LC_rtio(rep_ratio_counter))&&isfinite(group2_SA_rtio(rep_ratio_counter))...
            &&abs(group1_rep_zscore)<=1.65&&abs(group2_rep_zscore)<=1.65||isnan(group1_rep_zscore)&&isnan(group2_rep_zscore)&&group1_LC_1(zscore_counter)>setting.OD_threshuold&&group1_LC_2(zscore_counter)>setting.OD_threshuold
        
        
        ok_list(zscore_counter,1)=1;
    else
        ok_list(zscore_counter,1)=0;
    end
end
% next to convert the standart score to plate normalized valeus
plates_in_list=unique(plate_list);
% calculate the mean and standart deviation for all read in each plate the
% pass the oklist gate for both groups

if setting.use_madian==0
    for processing_counter=1:length(plates_in_list)
        if setting.use_abselute_deviation==1
            stdi_group1_LC_1(processing_counter)=mad(group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_LC_2(processing_counter)=mad(group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_SA_1(processing_counter)=mad(group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_SA_2(processing_counter)=mad(group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_LC_1(processing_counter)=mad(group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_LC_2(processing_counter)=mad(group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_SA_1(processing_counter)=mad(group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_SA_2(processing_counter)=mad(group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_LC(processing_counter)=mad(cat(1,group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group2_LC(processing_counter)=mad(cat(1,group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group1_SA(processing_counter)=mad(cat(1,group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group2_SA(processing_counter)=mad(cat(1,group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        else
            stdi_group1_LC_1(processing_counter)=nanstd(group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_LC_2(processing_counter)=nanstd(group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_SA_1(processing_counter)=nanstd(group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_SA_2(processing_counter)=nanstd(group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_LC_1(processing_counter)=nanstd(group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_LC_2(processing_counter)=nanstd(group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_SA_1(processing_counter)=nanstd(group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group2_SA_2(processing_counter)=nanstd(group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
            stdi_group1_LC(processing_counter)=nanstd(cat(1,group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group2_LC(processing_counter)=nanstd(cat(1,group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group1_SA(processing_counter)=nanstd(cat(1,group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
            stdi_group2_SA(processing_counter)=nanstd(cat(1,group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        end
        meani_group1_LC_1(processing_counter)=nanmean(group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_LC_2(processing_counter)=nanmean(group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_SA_1(processing_counter)=nanmean(group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_SA_2(processing_counter)=nanmean(group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_LC_1(processing_counter)=nanmean(group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_LC_2(processing_counter)=nanmean(group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_SA_1(processing_counter)=nanmean(group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_SA_2(processing_counter)=nanmean(group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_LC(processing_counter)=nanmean(cat(1,group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        meani_group2_LC(processing_counter)=nanmean(cat(1,group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        meani_group1_SA(processing_counter)=nanmean(cat(1,group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        meani_group2_SA(processing_counter)=nanmean(cat(1,group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
    end
elseif setting.use_madian==1
    for processing_counter=1:length(plates_in_list)
        meani_group1_LC_1(processing_counter)=nanmedian(group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_LC_2(processing_counter)=nanmedian(group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_SA_1(processing_counter)=nanmedian(group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group1_SA_2(processing_counter)=nanmedian(group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        stdi_group1_LC_1(processing_counter)=mad(group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group1_LC_2(processing_counter)=mad(group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group1_SA_1(processing_counter)=mad(group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group1_SA_2(processing_counter)=mad(group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        meani_group2_LC_1(processing_counter)=nanmedian(group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_LC_2(processing_counter)=nanmedian(group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_SA_1(processing_counter)=nanmedian(group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1));
        meani_group2_SA_2(processing_counter)=nanmedian(group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1));
        stdi_group2_LC_1(processing_counter)=mad(group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group2_LC_2(processing_counter)=mad(group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group2_SA_1(processing_counter)=mad(group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        stdi_group2_SA_2(processing_counter)=mad(group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1),1);
        meani_group1_LC(processing_counter)=nanmedian(cat(1,group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        stdi_group1_LC(processing_counter)=mad(cat(1,group1_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)),1);
        meani_group2_LC(processing_counter)=nanmedian(cat(1,group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        stdi_group2_LC(processing_counter)=mad(cat(1,group2_LC_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_LC_2(plate_list==plates_in_list(processing_counter)&ok_list==1)),1);
        meani_group1_SA(processing_counter)=nanmedian(cat(1,group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        stdi_group1_SA(processing_counter)=mad(cat(1,group1_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group1_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)),1);
        meani_group2_SA(processing_counter)=nanmedian(cat(1,group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)));
        stdi_group2_SA(processing_counter)=mad(cat(1,group2_SA_1(plate_list==plates_in_list(processing_counter)&ok_list==1),group2_SA_2(plate_list==plates_in_list(processing_counter)&ok_list==1)),1);
    end
end
% get the global mean and standart dev for all experament re scaling if
% selected

% maybe add aoption to rescale to only the plate valuses but for all musherments
if setting.use_madian==0
    if setting.use_abselute_deviation==0
        stdi_group1_LC_1_glob=nanstd(group1_LC_1(ok_list==1));
        stdi_group1_SA_1_glob=nanstd(group1_SA_1(ok_list==1));
        stdi_group2_LC_1_glob=nanstd(group2_LC_1(ok_list==1));
        stdi_group2_SA_1_glob=nanstd(group2_SA_1(ok_list==1));
    elseif setting.use_abselute_deviation==1
        stdi_group1_LC_1_glob=mad(group1_LC_1(ok_list==1));
        stdi_group1_SA_1_glob=mad(group1_SA_1(ok_list==1));
        stdi_group2_LC_1_glob=mad(group2_LC_1(ok_list==1));
        stdi_group2_SA_1_glob=mad(group2_SA_1(ok_list==1));
    end
    meani_group1_LC_1_glob=nanmean(group1_LC_1(ok_list==1));
    meani_group1_SA_1_glob=nanmean(group1_SA_1(ok_list==1));
    meani_group2_LC_1_glob=nanmean(group2_LC_1(ok_list==1));
    meani_group2_SA_1_glob=nanmean(group2_SA_1(ok_list==1));
elseif setting.use_madian==1
    meani_group1_LC_1_glob=nanmedian(group1_LC_1(ok_list==1));
    meani_group1_SA_1_glob=nanmedian(group1_SA_1(ok_list==1));
    stdi_group1_LC_1_glob=mad(group1_LC_1(ok_list==1),1);
    stdi_group1_SA_1_glob=mad(group1_SA_1(ok_list==1),1);
    meani_group2_LC_1_glob=nanmedian(group2_LC_1(ok_list==1));
    meani_group2_SA_1_glob=nanmedian(group2_SA_1(ok_list==1));
    stdi_group2_LC_1_glob=mad(group2_LC_1(ok_list==1),1);
    stdi_group2_SA_1_glob=mad(group2_SA_1(ok_list==1),1);
end

if setting.LC_normalization_factor==1
    LC_normalization_factor_value=1000;
elseif setting.LC_normalization_factor==0
    LC_normalization_factor_value=(meani_group1_LC_1_glob+meani_group2_LC_1_glob)/2;
end
%re beild the out table with normalized number

temp_processed_table=cell(size(data_no_beckgrund,1),15);
for final_counter1=1:size(data_no_beckgrund,1)
    final_counter2=final_counter1;
    if setting.LC_global_normalization==1
        if ok_list(final_counter1)==1
            %get the pate number and cordination from sumtable
            temp_processed_table{final_counter2,1}=data_no_beckgrund{final_counter2,1};
            temp_processed_table{final_counter2,2}=data_no_beckgrund{final_counter2,2};
            temp_processed_table{final_counter2,3}=setting.libcord.Gene(setting.libcord.cord384==data_no_beckgrund{final_counter2,2});
            temp_processed_table{final_counter2,4}=setting.libcord.ORF(setting.libcord.cord384==data_no_beckgrund{final_counter2,2});
            % identefay the plate number
            temp_plate=find(plates_in_list==plate_list(final_counter2));
            % calculate the standart score using the plate mean and stds and re
            % scale to the global mean and stds for the group for each repet
            % and each muserment and avarage the result
            % normalize the SA read using the LC
            if setting.use_OD_read ==1
            temp_processed_table{final_counter2,5}=data_no_beckgrund{final_counter2,3};
            
            temp_processed_table{final_counter2,9}=data_no_beckgrund{final_counter2,9};
            
            elseif setting.use_OD_read ==0
            temp_processed_table{final_counter2,5}=(((((data_no_beckgrund{final_counter2,3}-meani_group1_LC_1(temp_plate))/stdi_group1_LC_1(temp_plate))*stdi_group1_LC_1_glob)+meani_group1_LC_1_glob)+...
                (((data_no_beckgrund{final_counter2,5}-meani_group1_LC_2(temp_plate))/stdi_group1_LC_2(temp_plate))*stdi_group1_LC_1_glob)+meani_group1_LC_1_glob)/2;
            
            temp_processed_table{final_counter2,9}=(((((data_no_beckgrund{final_counter2,9}-meani_group2_LC_1(temp_plate))/stdi_group2_LC_1(temp_plate))*stdi_group2_LC_1_glob)+meani_group2_LC_1_glob)+...
                (((data_no_beckgrund{final_counter2,11}-meani_group2_LC_2(temp_plate))/stdi_group2_LC_2(temp_plate))*stdi_group2_LC_1_glob)+meani_group2_LC_1_glob)/2;
            end
            temp_processed_table{final_counter2,6}=(((((data_no_beckgrund{final_counter2,6}-meani_group1_SA_2(temp_plate))/stdi_group1_SA_2(temp_plate))*stdi_group1_SA_1_glob)+meani_group1_SA_1_glob)+...
                (((data_no_beckgrund{final_counter2,4}-meani_group1_SA_1(temp_plate))/stdi_group1_SA_1(temp_plate))*stdi_group1_SA_1_glob)+meani_group1_SA_1_glob)/2;
            
            temp_processed_table{final_counter2,7}=temp_processed_table{final_counter2,6}*(LC_normalization_factor_value/temp_processed_table{final_counter2,5});
            
              
            temp_processed_table{final_counter2,10}=(((((data_no_beckgrund{final_counter2,10}-meani_group2_SA_1(temp_plate))/stdi_group2_SA_1(temp_plate))*stdi_group2_SA_1_glob)+meani_group2_SA_1_glob)+...
                (((data_no_beckgrund{final_counter2,12}-meani_group2_SA_2(temp_plate))/stdi_group2_SA_2(temp_plate))*stdi_group2_SA_1_glob)+meani_group2_SA_1_glob)/2;
            
            temp_processed_table{final_counter2,11}=temp_processed_table{final_counter2,10}*(LC_normalization_factor_value/temp_processed_table{final_counter2,9});
            % remove nagative noms and calculate the group1 to group2 ratio
            if temp_processed_table{final_counter2,7}<0||temp_processed_table{final_counter2,11}<0
                temp_processed_table{final_counter2,13}=0;
            else
                temp_processed_table{final_counter2,13}=(temp_processed_table{final_counter2,7}/temp_processed_table{final_counter2,11});
            end
            %if the ratio is higer the the treshold(th2) set the group preferanse if not set it to same
            if temp_processed_table{final_counter2,6}<setting.signal_threshuold||temp_processed_table{final_counter2,10}<setting.signal_threshuold
                temp_processed_table{final_counter2,14}='below TH';
            elseif  temp_processed_table{final_counter2,13}>setting.ratio_treshuld
                temp_processed_table{final_counter2,14}=group1;
            elseif temp_processed_table{final_counter2,13}<(1/setting.ratio_treshuld)
                temp_processed_table{final_counter2,14}=group2;
            else
                temp_processed_table{final_counter2,14}='same';
            end
        end
    elseif setting.LC_global_normalization==0
        if ok_list(final_counter1)==1
            %get the pate number and cordination from sumtable
            temp_processed_table{final_counter2,1}=data_no_beckgrund{final_counter2,1};
            temp_processed_table{final_counter2,2}=data_no_beckgrund{final_counter2,2};
            temp_processed_table{final_counter2,3}=setting.libcord.Gene(setting.libcord.cord384==data_no_beckgrund{final_counter2,2});
            temp_processed_table{final_counter2,4}=setting.libcord.ORF(setting.libcord.cord384==data_no_beckgrund{final_counter2,2});
            
            % identefay the plate number
            temp_plate=find(plates_in_list==plate_list(final_counter2));
            % calculate the standart score using the plate mean and stds and re
            % scale to the global mean and stds for the group for each repet
            % and each muserment and avarage the result
            % normalize the SA read using the LC
            if setting.use_OD_read ==1
            temp_processed_table{final_counter2,5}=data_no_beckgrund{final_counter2,3};
            
            temp_processed_table{final_counter2,9}=data_no_beckgrund{final_counter2,9};
            
            elseif setting.use_OD_read ==0
                temp_processed_table{final_counter2,5}=(((((data_no_beckgrund{final_counter2,3}-meani_group1_LC_1(temp_plate))/stdi_group1_LC_1(temp_plate))*stdi_group1_LC_1(temp_plate))+meani_group1_LC_1(temp_plate))+...
                    (((data_no_beckgrund{final_counter2,5}-meani_group1_LC_2(temp_plate))/stdi_group1_LC_2(temp_plate))*stdi_group1_LC_1(temp_plate))+meani_group1_LC_1(temp_plate))/2;
                
                temp_processed_table{final_counter2,9}=(((((data_no_beckgrund{final_counter2,9}-meani_group2_LC_1(temp_plate))/stdi_group2_LC_1(temp_plate))*stdi_group2_LC_1(temp_plate))+meani_group2_LC_1(temp_plate))+...
                    (((data_no_beckgrund{final_counter2,11}-meani_group2_LC_2(temp_plate))/stdi_group2_LC_2(temp_plate))*stdi_group2_LC_1(temp_plate))+meani_group2_LC_1(temp_plate))/2;
            end
            
                        
            temp_processed_table{final_counter2,6}=(((((data_no_beckgrund{final_counter2,6}-meani_group1_SA_2(temp_plate))/stdi_group1_SA_2(temp_plate))*stdi_group1_SA_1(temp_plate))+meani_group1_SA_1(temp_plate))+...
                (((data_no_beckgrund{final_counter2,4}-meani_group1_SA_1(temp_plate))/stdi_group1_SA_1(temp_plate))*stdi_group1_SA_1(temp_plate))+meani_group1_SA_1(temp_plate))/2;
            
            
            temp_processed_table{final_counter2,7}=temp_processed_table{final_counter2,6}*(LC_normalization_factor_value/temp_processed_table{final_counter2,5});
            
            
            temp_processed_table{final_counter2,10}=(((((data_no_beckgrund{final_counter2,10}-meani_group2_SA_1(temp_plate))/stdi_group2_SA_1(temp_plate))*stdi_group2_SA_1(temp_plate))+meani_group2_SA_1(temp_plate))+...
                (((data_no_beckgrund{final_counter2,12}-meani_group2_SA_2(temp_plate))/stdi_group2_SA_2(temp_plate))*stdi_group2_SA_1(temp_plate))+meani_group2_SA_1(temp_plate))/2;
            temp_processed_table{final_counter2,11}=temp_processed_table{final_counter2,10}*(LC_normalization_factor_value/temp_processed_table{final_counter2,9});
            
            % remove nagative noms and calculate the group1 to group2 ratio
            if temp_processed_table{final_counter2,7}<0||temp_processed_table{final_counter2,11}<0
                temp_processed_table{final_counter2,13}=0;
            else
                temp_processed_table{final_counter2,13}=(temp_processed_table{final_counter2,7}/temp_processed_table{final_counter2,11});
            end
            %if the ratio is higer the the treshold(th2) set the group preferanse if not set it to same
            if temp_processed_table{final_counter2,6}<setting.signal_threshuold||temp_processed_table{final_counter2,10}<setting.signal_threshuold
                temp_processed_table{final_counter2,14}='below TH';
            elseif  temp_processed_table{final_counter2,13}>setting.ratio_treshuld
                temp_processed_table{final_counter2,14}=group1;
            elseif temp_processed_table{final_counter2,11}<(1/setting.ratio_treshuld)
                temp_processed_table{final_counter2,14}=group2;
                
            else
                temp_processed_table{final_counter2,14}='same';
            end
        end
    end
end
for temp_table_counter=2:size(temp_processed_table,1)
    if ~(isempty(temp_processed_table{temp_table_counter,5})||temp_processed_table{temp_table_counter,5}<0||isempty(temp_processed_table{temp_table_counter,6})||temp_processed_table{temp_table_counter,6}<0||isempty(temp_processed_table{temp_table_counter,7})||temp_processed_table{temp_table_counter,7}<0||isempty(temp_processed_table{temp_table_counter,9})||temp_processed_table{temp_table_counter,9}<0||isempty(temp_processed_table{temp_table_counter,10})||temp_processed_table{temp_table_counter,10}<0||isempty(temp_processed_table{temp_table_counter,11})||temp_processed_table{temp_table_counter,11}<0||temp_processed_table{temp_table_counter,3}=="None")
        processed_data=vertcat(processed_data,{temp_processed_table{temp_table_counter,1:14}});
    end
end

for noncounter=1:length(data_no_beckgrund)
    temp_sum_table1{noncounter+1,1}=data_no_beckgrund{noncounter,1};
    temp_sum_table1{noncounter+1,2}=data_no_beckgrund{noncounter,2};
    temp_sum_table1{noncounter+1,3}=setting.libcord.Gene(setting.libcord.cord384==data_no_beckgrund{noncounter,2});
    temp_sum_table1{noncounter+1,4}=setting.libcord.ORF(setting.libcord.cord384==data_no_beckgrund{noncounter,2});
end


data_no_beckgrundff=cat(1,nonprocessed_data,data_no_beckgrund);
nonprocessed_table_out=[temp_sum_table1(:,1:4),data_no_beckgrundff(:,3:end)];

disp('saveing the output files')
if setting.save_tabels==1
    
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'processed_data'))
    if  setting.save_nonnormalaized_data==1
        nonprocesd_file_name=strcat(folder_path,filesep,'analesis_files',filesep,'processed_data',filesep,setting.out_file_initial,'_nonnormalaized_data_',out_name_string,'.csv');
        writecell(nonprocessed_table_out,nonprocesd_file_name)
    end
    
    procesd_file_name=strcat(folder_path,filesep,'analesis_files',filesep,'processed_data',filesep,setting.out_file_initial,'_procesd_data_',out_name_string,'.csv');
    writecell(processed_data,procesd_file_name)
end

end