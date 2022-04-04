if exist('settinghash','var') == 0
    settinghash=[];
end
[setting,settinghash] = settingdef(settinghash);

if exist('row_data','var') == 0 && exist('folder_path','var') == 0
    disp('extrecting data from images')
    [row_data,folder_path] = extract_data_from_images_and_ROI(setting);
else
    disp('using preload row data')
end

if exist('processed_data','var') == 0 && exist('nonprocessed_table_out','var')== 0&& exist('group1','var') == 0 && exist('group2','var')== 0&& exist('plates_in_list','var')== 0&& exist('out_name_string','var')== 0||setting.test_mode==1
[processed_data,nonprocessed_table_out,group1,group2]=...
    dotblot_analisis(setting,row_data,folder_path);
else
    disp('using preload analisis data')
end
disp('ploting the results')
if setting.plot==1
    plotall(processed_data,nonprocessed_table_out,setting.plot_sum,setting.plot_signal_disterbution,setting.plot_corolation,setting.plot_ratio,setting.save_figs,folder_path)
end



disp('sequnce analesis')
if setting.sequnce_analesis==1
   [intrestingaa,aalist]=sequnce_analesis(processed_data,group1,group2,folder_path,setting);
end