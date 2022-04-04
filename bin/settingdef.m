function [setting,settinghash] = settingdef(settinghashold)
setting.show_figs=0;
setting.save_tabels=1;
setting.save_nonnormalaized_data=1;
setting.plot=0;
setting.save_figs=1;
setting.sequnce_analesis=1;
setting.save_sequance_analesis_data=1;
setting.do_hlixwhell=0;
setting.do_seqlogo=0;

%general settings
setting.install_path=fileparts(fileparts(mfilename('fullpath')));
setting.lib_path=strcat(setting.install_path,filesep,'lib');
setting.libcord = getlibcord(strcat(setting.lib_path,filesep,'libcord.csv'));
setting.topologyandseq = gettopology(strcat(setting.lib_path,filesep,'topologyandseq.csv'));
%setting.SGD_DB = SGD_DB_import(strcat(setting.lib_path,filesep,'SGDDB.csv'));


%difine tresholds and normalization type
setting.pre_normalization_type=0;%0= no prenormalization, 1= minmax scailing 
setting.LC_normalization_factor=1;%lc re scaling factor 0=1k factor 1= global LC mean 2(toadd)=global LC median
setting.LC_global_normalization=1;%0=per plate normalization 1= global plate normalization
setting.use_madian=0;%0=use mean 1=use meadian
setting.use_abselute_deviation=0;%0=use standartdeviation 1=use abeslut deviation (only works with median statistics)
setting.use_dot_stat_median=1;%statistic to extract from the image 0=mean,1=median
setting.use_OD_read=1;
setting.min_signal_threshuold=0;% valeus below this int will be excloded
setting.signal_threshuold=1000;
setting.ratio_treshuld=2;% retio factor for selecting the hits
setting.OD_threshuold=0.5;
%ploting
setting.plot_sum=0;
setting.plot_signal_disterbution=0;
setting.plot_corolation=0;
setting.plot_ratio=0;

%tabels
setting.out_file_initial='outfile';


%sequnce analesis

setting.get_TMD=0;
setting.get_SS=1;
setting.get_cyto=0;
setting.get_luminal=0;
setting.min_TMD_size=3;
setting.max_TMD_size=25;
setting.TMD_to_check=1;%'maxim';
setting.set_TMD_size=20;%0= get max length from predection
setting.set_TMD_size_exect=1;%1=exect set_TMD_size 0=more if exist
setting.hits_uni=0;% 1 get uniun,0 get intersect
setting.helixwheel_sorting=1;%  sort the pie according to the propertis scale 0= so sorting presentetion acording to the abc
setting.helixwheel_selected_fun=5;%the nomber of the function !! nned to change to name !!
setting.helix_number=3;
setting.save_intermidiet_whells=0;
setting.save_intermidiet_seqlogo=0;
%sequnce analesis functions to run
setting.get_functions_mean=0;
setting.spesific_AA=0;
setting.func = {...

%stractural propertis related
@aaaverageburied;
@aaaverageflexibility;
@aabulkiness;
@aaburiedresidues;

%chimical propertis related
@aahphob_argos;
@aahphob_black;
@aahphob_breese;
@aahphob_chothia;

@aahphob_doolittle;
@aahphob_eisenberg;
@aahphob_fauchere;
@aahphob_guy;
@aahphob_janin;
@aahphob_leo;
@aahphob_manavalan;
@aahphob_miyazawa;
@aahphob_mobility;
@aahphob_parker;
@aahphob_rose;

@aahphob_roseman;
@aahphob_sweet;
@aahphob_welling;
@aahphob_wilson;
@aahphob_wolfenden;
@aahphob_woods;

@aapolaritygrantham;
@aapolarityzimmerman;

@AAIndex

@aaratioside;
@aarecognitionfactors;
@aarefractivity;

@aahphob_ph3_4;
@aahphob_ph7_5;




%not intresting

%alpha_helix related
@aaalpha_helixfasman;
@aaalpha_helixlevitt;
@aaalpha_helixroux;

%coil related
@aacoilroux;

%beta_sheet related
@aaantiparallelbeta_strand;
@aabeta_sheetfasman;
@aabeta_sheetlevitt;
@aabeta_sheetroux;
@aabeta_turnfasman;
@aabeta_turnlevitt;
@aabeta_turnroux;
@aaparallelbeta_strand;
@aatotalbeta_strand;

@aarelativemutability;
@aamolecularweight;

@aahplc2_1;
@aahplc7_4;
@aahplchfba;
@aahplctfa;

};
setting.aa = ['A'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'K'; 'L'; 'M'; 'N'; 'P'; 'Q'; 'R'; 'S'; 'T'; 'V'; 'W'; 'Y'];

%checing if statistics recalcolation is requierd
setting.test_mode=0;
settinghash=DataHash(rmfield(setting,'test_mode'));

if isempty(settinghashold)
    disp('settings struct ctreated')
elseif ~strcmp(settinghash,settinghashold)
    setting.test_mode=1;
    disp('settings changed recalculating statistics')
    
else
    setting.test_mode=0;
    disp('settings dident changed')

end
end
