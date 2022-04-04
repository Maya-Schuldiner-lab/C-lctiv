function [row_data,folder_path] = extract_data_from_images_and_ROI(setting)
%read librery cordination and list images
warning('off', 'MATLAB:MKDIR:DirectoryExists');
folder_path = uigetdir;
filesAndFolders=dir(fullfile(folder_path, '*.tif'));
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));
odreads=readtable(strcat(folder_path,filesep,'odreads.csv'));
%% make well
%prelocation

numtocord=cell([384 2]);
abclist='ABCDEFGHIJKLMNOP';
well_counter=1;
for Col=1:length(abclist)
    for Row=1:24
        numtocord{well_counter,1}=well_counter;
        numtocord{well_counter,2}=strcat(abclist(Col),string(Row));
        well_counter=well_counter+1;
    end
end

%%
header_part_1{1,1}='plate_number';
header_part_1{1,2}='sample_cord';
header_part_1{1,3}='ORF';
header_part_1{1,4}='Gene_name';
%muserments

%%
%data extraction
%prelocation

file_name_and_ext=cell(1,length(filesInDir));
ROI_file_name=cell(1,length(filesInDir));
ROI_data=cell(1,length(filesInDir));
file_name_vector=cell(1,length(filesInDir));
sample_name=cell(1,length(filesInDir));
sample_plate_number=cell(1,length(filesInDir));
ch1=cell(1,length(filesInDir));
ch2=cell(1,length(filesInDir));
groups=cell(1,length(filesInDir));
singel_dot_ch2=cell(1,384);
singel_dot_ch1=cell(1,384);
BW=cell(1,384);
maskedImage=cell(1,384);
BW_inv=cell(1,384);
maskedImage_inv=cell(1,384);

textwaitbar(0, length(filesInDir), 'reading the dots')
for row_data_files=1:length(filesInDir)
    
    file_name_and_ext{row_data_files}=split(filesInDir(row_data_files).name,'.');
    ROI_file_name{row_data_files}=strcat(filesInDir(row_data_files).folder,filesep,string(file_name_and_ext{row_data_files}(1)),'_Results Documents mean',filesep,string(file_name_and_ext{row_data_files}(1)),'_roi.zip');
    ROI_data{row_data_files}=ReadImageJROI(ROI_file_name{row_data_files});
    file_name_vector{row_data_files}=split(filesInDir(row_data_files).name,' ');
    sample_name{row_data_files}=split(file_name_vector{row_data_files}{2},'.');
    
    plate1536_number=sample_name{row_data_files}{1}(2);
    plate1536_quad=str2double(sample_name{row_data_files}{1}(4));
    Repeat(row_data_files)=str2double(sample_name{row_data_files}{1}(6));
    
    % convert the 1536 plate number and quad number to 384 plat number (in
    % v.01 is build for avitag and bira laib so max 1536 number is 4
    if plate1536_number=='1'
        sample_plate_number{row_data_files}=plate1536_quad;
    elseif plate1536_number=='2'
        sample_plate_number{row_data_files}=4+plate1536_quad;
    elseif plate1536_number=='3'
        sample_plate_number{row_data_files}=8+plate1536_quad;
    elseif plate1536_number=='4'
        sample_plate_number{row_data_files}=12+plate1536_quad;
    end
    %read images ch1 is 680nm and ch2 is 800nm
    ch1{row_data_files}=imread([filesInDir(row_data_files).folder,filesep,filesInDir(row_data_files).name],1);
    ch2{row_data_files}=imread([filesInDir(row_data_files).folder,filesep,filesInDir(row_data_files).name],2);
    
    
    group_id=string(file_name_vector{row_data_files}(1));
    groups{row_data_files}=group_id;
    plate_starting_counter=(sample_plate_number{row_data_files}*384)+2-384;
    Repeat_index=(Repeat(row_data_files)-1)*41;
    for well_counter=1:384
        %rectangle('Position',[roidata{1,row_data_files}{1,i}.vnRectBounds(2)-3 roidata{1,row_data_files}{1,i}.vnRectBounds(1)-3 22 22],'EdgeColor','w')
        
        singel_dot_ch1{well_counter}=imcrop(ch1{row_data_files},[ROI_data{1,row_data_files}{1,well_counter}.vnRectBounds(2)-3 ROI_data{1,row_data_files}{1,well_counter}.vnRectBounds(1)-3 22 22]);
        singel_dot_ch2{well_counter}=imcrop(ch2{row_data_files},[ROI_data{1,row_data_files}{1,well_counter}.vnRectBounds(2)-3 ROI_data{1,row_data_files}{1,well_counter}.vnRectBounds(1)-3 22 22]);
        [BW{well_counter},maskedImage{well_counter},BW_inv{well_counter},maskedImage_inv{well_counter}] = segmentdot(singel_dot_ch2{well_counter});
        
        %metadata
        semple_data.(group_id){plate_starting_counter,1}=string(sample_plate_number{row_data_files});
        semple_data.(group_id){plate_starting_counter,2}=strcat(string(sample_plate_number{row_data_files}),numtocord{well_counter,2});
        semple_data.(group_id){plate_starting_counter,3}=setting.libcord.ORF(setting.libcord.cord384==strcat(string(sample_plate_number{row_data_files}),numtocord{well_counter,2}));
        semple_data.(group_id){plate_starting_counter,4}=setting.libcord.Gene(setting.libcord.cord384==strcat(string(sample_plate_number{row_data_files}),numtocord{well_counter,2}));
        
        %muserments
        dch1f=double(singel_dot_ch1{well_counter}(BW{well_counter}));
        measurement_data.(group_id){plate_starting_counter,1+Repeat_index}=std(dch1f);
        measurement_data.(group_id){plate_starting_counter,2+Repeat_index}=mean(dch1f);
        measurement_data.(group_id){plate_starting_counter,3+Repeat_index}=median(dch1f);
        measurement_data.(group_id){plate_starting_counter,4+Repeat_index}=min(dch1f);
        measurement_data.(group_id){plate_starting_counter,5+Repeat_index}=max(dch1f);

       
        dch1b=double(singel_dot_ch1{well_counter}(BW_inv{well_counter}));
        measurement_data.(group_id){plate_starting_counter,6+Repeat_index}=std(dch1b);
        measurement_data.(group_id){plate_starting_counter,7+Repeat_index}=mean(dch1b);
        measurement_data.(group_id){plate_starting_counter,8+Repeat_index}=median(dch1b);
        measurement_data.(group_id){plate_starting_counter,9+Repeat_index}=min(dch1b);
        measurement_data.(group_id){plate_starting_counter,10+Repeat_index}=max(dch1b);
        
        
        
        dch2f=double(singel_dot_ch2{well_counter}(BW{well_counter}));
        measurement_data.(group_id){plate_starting_counter,11+Repeat_index}=std(dch2f);
        measurement_data.(group_id){plate_starting_counter,12+Repeat_index}=mean(dch2f);
        measurement_data.(group_id){plate_starting_counter,13+Repeat_index}=median(dch2f);
        measurement_data.(group_id){plate_starting_counter,14+Repeat_index}=min(dch2f);
        measurement_data.(group_id){plate_starting_counter,15+Repeat_index}=max(dch2f);
        
        dch2b=double(singel_dot_ch2{well_counter}(BW_inv{well_counter}));
        measurement_data.(group_id){plate_starting_counter,16+Repeat_index}=std(dch2b);
        measurement_data.(group_id){plate_starting_counter,17+Repeat_index}=mean(dch2b);
        measurement_data.(group_id){plate_starting_counter,18+Repeat_index}=median(dch2b);
        measurement_data.(group_id){plate_starting_counter,19+Repeat_index}=min(dch2b);
        measurement_data.(group_id){plate_starting_counter,20+Repeat_index}=max(dch2b);
       %% > b680 + 1        
        measurement_data.(group_id){plate_starting_counter,21+Repeat_index}=sum(dch1f>(mean(dch1b)+std(dch1b)))*100/sum(sum(BW_inv{well_counter}));
        %% > b680 + 2
        measurement_data.(group_id){plate_starting_counter,22+Repeat_index}=sum(dch1f>(mean(dch1b)+(2*std(dch1b))))*100/sum(sum(BW_inv{well_counter}));
        
        %% > b800 + 1
        measurement_data.(group_id){plate_starting_counter,23+Repeat_index}=sum(dch2f>(mean(dch2b)+std(dch2b)))*100/sum(sum(BW{well_counter}));
        %% > b800 + 2
        measurement_data.(group_id){plate_starting_counter,24+Repeat_index}=sum(dch2f>(mean(dch2b)+(2*std(dch2b))))*100/sum(sum(BW{well_counter}));
        %% ratio of mean 
        measurement_data.(group_id){plate_starting_counter,25+Repeat_index}= mean((dch1f-mean(dch1b)))/mean((dch2f-mean(dch2b)));
        %% ratio of median 
        measurement_data.(group_id){plate_starting_counter,26+Repeat_index}= median((dch1f-mean(dch1b)))/median((dch2f-mean(dch2b)));
        %% Log Ratio
        measurement_data.(group_id){plate_starting_counter,27+Repeat_index}=log2(median((dch1f-mean(dch1b)))/median((dch2f-mean(dch2b))));
        %% Mean of Ratios
        measurement_data.(group_id){plate_starting_counter,28+Repeat_index}= mean((dch2f-mean(dch2b))./(dch1f-mean(dch1b)));
        %% Median of Ratios
        measurement_data.(group_id){plate_starting_counter,29+Repeat_index}=median((dch2f-mean(dch2b))./(dch1f-mean(dch1b)));
        %% Ratios SD
        measurement_data.(group_id){plate_starting_counter,30+Repeat_index}=std((dch2f-mean(dch2b))./(dch1f-mean(dch1b)));
        %% F Pixels 
        measurement_data.(group_id){plate_starting_counter,31+Repeat_index}=sum(sum(BW{well_counter}));
        %% B Pixels 
        measurement_data.(group_id){plate_starting_counter,32+Repeat_index}=sum(sum(BW_inv{well_counter}));
        %% F1 Median - B1 
        measurement_data.(group_id){plate_starting_counter,33+Repeat_index}=median(dch1f)-median(dch1b);
        %% F2 Median - B2
        measurement_data.(group_id){plate_starting_counter,34+Repeat_index}=median(dch2f)-median(dch2b);
         %% F1 Mean - B1
        measurement_data.(group_id){plate_starting_counter,35+Repeat_index}=mean(dch1f)-median(dch1b);

        %% F2 Mean - B2
        measurement_data.(group_id){plate_starting_counter,36+Repeat_index}=mean(dch2f)-median(dch2b);
        %% SNR 1
        measurement_data.(group_id){plate_starting_counter,37+Repeat_index}=(mean(dch1f)-mean(dch1b))/std(dch1b);
        %% SNR 2
        measurement_data.(group_id){plate_starting_counter,38+Repeat_index}=(mean(dch2f)-mean(dch2b))/std(dch2b);
        %% F1 Total Intensity
        measurement_data.(group_id){plate_starting_counter,39+Repeat_index}=sum(dch1f);
        %% F2 Total Intensity
        measurement_data.(group_id){plate_starting_counter,40+Repeat_index}=sum(dch2f);
        measurement_data.(group_id){plate_starting_counter,41+Repeat_index}=odreads.(group_id)(odreads.cord==strcat(string(sample_plate_number{row_data_files}),numtocord{well_counter,2}));
        plate_starting_counter=plate_starting_counter+1;
        
    end
    
    %file counter
    textwaitbar(row_data_files, length(filesInDir), 'reading the dots')
    
end
%build headers for the reapet number
for dd=string(unique(Repeat))
    header_template{1,1}=strcat('F680_SD_R',dd);
    header_template{1,2}=strcat('F680_mean_R',dd);
    header_template{1,3}=strcat('F680_median_R',dd);
    header_template{1,4}=strcat('F680_min_R',dd);
    header_template{1,5}=strcat('F680_max_R',dd);
    header_template{1,6}=strcat('B680_SD_R',dd);
    header_template{1,7}=strcat('B680_mean_R',dd);
    header_template{1,8}=strcat('B680_median_R',dd);
    header_template{1,9}=strcat('B680_min_R',dd);
    header_template{1,10}=strcat('B680_max_R',dd);
    header_template{1,11}=strcat('F800_SD_R',dd);
    header_template{1,12}=strcat('F800_mean_R',dd);
    header_template{1,13}=strcat('F800_median_R',dd);
    header_template{1,14}=strcat('F800_min_R',dd);
    header_template{1,15}=strcat('F800_max_R',dd);
    header_template{1,16}=strcat('B800_SD_R',dd);
    header_template{1,17}=strcat('B800_mean_R',dd);
    header_template{1,18}=strcat('B800_median_R',dd);
    header_template{1,19}=strcat('B800_min_R',dd);
    header_template{1,20}=strcat('B800_max_R',dd);
    header_template{1,21}=strcat('f680_more_b680_plas_1_R',dd);
    header_template{1,22}=strcat('f680_more_b680_plas_2_R',dd);        
    header_template{1,23}=strcat('f800_more_b800_plas_1_R',dd);
    header_template{1,24}=strcat('f800_more_b800_plas_2_R',dd);
    header_template{1,25}=strcat('ratio_of_mean_R',dd); 
    header_template{1,26}=strcat('ratio_of_median_R',dd); 
    header_template{1,27}=strcat('Log_Ratio_R',dd);
    header_template{1,28}=strcat('Mean_of_Ratios_R',dd);
    header_template{1,29}=strcat('Median_of_Ratios_R',dd);
    header_template{1,30}=strcat('Ratios_SD_R',dd);
    header_template{1,31}=strcat('F_Pixels_number_R',dd);
    header_template{1,32}=strcat('B_Pixels_number_R',dd);
    header_template{1,33}=strcat('F680_Median_minus_Median_R',dd); 
    header_template{1,34}=strcat('F800_Median_minus_Median_R',dd);
    header_template{1,35}=strcat('F680_Mean_minus_Median_R',dd);
    header_template{1,36}=strcat('F800_Mean_minus_Median_R',dd);
    header_template{1,37}=strcat('SNR_680_R',dd);
    header_template{1,38}=strcat('SNR_800_R',dd);
    header_template{1,39}=strcat('F680_Total_Intensity_R',dd);
    header_template{1,40}=strcat('F800_Total_Intensity_R',dd);
    header_template{1,41}=strcat('OD_Read_R',dd);

    header_part_1=cat(2,header_part_1,header_template);
end


groups_to_save=unique(string(groups));
for group_index2=1:length(groups_to_save)
    all_data.(groups_to_save(group_index2))=cat(2,semple_data.(groups_to_save(group_index2)),measurement_data.(groups_to_save(group_index2)));
    all_data_with_headers.(groups_to_save(group_index2))=cat(1,header_part_1,all_data.(groups_to_save(group_index2))(2:end,:));
    try
        mkdir(strcat(folder_path,filesep,'analesis_files'))
    catch
    end
    try
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'row_image_data'))
    catch
    end
    procesd_file_name=strcat(folder_path,filesep,'analesis_files',filesep,'row_image_data',filesep,'procesd_data_',groups_to_save(group_index2),'.csv');
    fprintf('saving row results for %s\n', string(groups_to_save(group_index2)))
    writecell(all_data_with_headers.(groups_to_save(group_index2)),procesd_file_name)
    row_data.(groups_to_save(group_index2))=cell2struct(all_data.(groups_to_save(group_index2))(2:end,:),[header_part_1{:}],2);
    %row_data.(groups_to_save(group_index2))=cell2table(all_data_with_headers.(groups_to_save(group_index2)));
end
ver_file_name=strcat(folder_path,filesep,'analesis_files',filesep,'row_image_data',filesep,'procesd_data_ver.mat');
save(ver_file_name,'row_data','folder_path')

end
