function [intrestingaa,aalist]=sequnce_analesis(processed_data,group1,group2,folder_path,setting)
disp('runing sequnce analysis')
%for internal pip data

hits{1,1}='ORF';
hits{1,2}='hit ratio';
hits{1,3}='hit';
chhc=2;
for wq=2:length(processed_data)
    if strcmp(processed_data{wq,14},group1)||strcmp(processed_data{wq,14},group2)||strcmp(processed_data{wq,14},'same')||strcmp(processed_data{wq,14},'all')
        hits{chhc,1}=processed_data{wq,4};
        hits{chhc,2}=processed_data{wq,13};
        hits{chhc,3}=processed_data{wq,14};
        chhc=chhc+1;
    end
end
%{
hits{1,1}='ORF';
hits{1,2}='hit type';
hits{1,3}='hit';
chhc=2;
for wq=2:length(xtry)
    if strcmp(xtry{wq,3},group1)||strcmp(xtry{wq,3},group2)||strcmp(xtry{wq,3},'same')||strcmp(xtry{wq,3},'all')
        hits{chhc,1}=xtry{wq,1};
        hits{chhc,2}=xtry{wq,2};
        hits{chhc,3}=xtry{wq,3};
        chhc=chhc+1;
    end
end
%}


allfilednames=setting.topologyandseq.Properties.VariableNames;
predictors=allfilednames(4:end)';
topologyandseqq = array2table(zeros(0,0));
warning('off','all')
for ct=2:length(hits)
    if ~isempty(find(setting.topologyandseq.ORF==upper(hits{ct,1}), 1))
        for fif=1:length(allfilednames)
            
            topologyandseqq.(allfilednames{fif})(ct-1)=setting.topologyandseq.(allfilednames{fif})(setting.topologyandseq.ORF==upper(hits{ct,1}));
        end
        topologyandseqq.hitname(ct-1)=string(hits{ct,3});
        topologyandseqq.hitratio(ct-1)=hits{ct,2};
    end
    textwaitbar(ct, length(hits), 'geting topology data')
    
end
warning('on','all')
warning('off', 'MATLAB:MKDIR:DirectoryExists');

for count=1:height(topologyandseqq)
    topologyandseqlog(count).seq =topologyandseqq.SequencefromSGD(count);
    topologyandseqlog(count).name =topologyandseqq.genename(count);
    topologyandseqlog(count).hitratio=topologyandseqq.hitratio(count);
    topologyandseqlog(count).hitname=topologyandseqq.hitname(count);

    for p=1:length(predictors)
        if setting.get_functions_mean==1
            topologyandseqlog(count).(predictors{p})=topologyandseqq.(predictors{p})(count);
            
        else
        try
            string1=topologyandseqq.(predictors{p})(count);
            if setting.get_luminal==1
            string1 = strrep(string1,'o','1');
            else
            string1 = strrep(string1,'o','0');
            end
            
            if setting.get_cyto==1
            string1 = strrep(string1,'i','1');
            else
            string1 = strrep(string1,'i','0');
            end
            
            if setting.get_TMD==1
                string1 = strrep(string1,'M','1');
            else
                string1 = strrep(string1,'M','0');
            end
            
            if setting.get_SS==1
                string1 = strrep(string1,'S','1');
            else
                string1 = strrep(string1,'S','0');
            end
            
            logstring=logical(char(string1(:))-'0');
            
            [labeledRegions, numberOfRegions] = bwlabel(bwareafilt(logstring,[setting.min_TMD_size,setting.max_TMD_size]));
            topologyandseqlog(count).TMDnum=numberOfRegions;
            if numberOfRegions>0
                %for tmdn=1:numberOfRegions
                
                tempseq=char(topologyandseqlog(count).seq);
                if setting.TMD_to_check=='maxim'
                idx= minmax(find(labeledRegions==numberOfRegions));
                else
                idx= minmax(find(labeledRegions==setting.TMD_to_check));
                end
                if idx(1)==1
                    idx(1)=2;
                end
                tmd_size=idx(2)-idx(1)+1;
                if setting.set_TMD_size==0
                    
                    topologyandseqlog(count).(predictors{p})=tempseq(idx(1):idx(2));
                elseif setting.set_TMD_size<=tmd_size||setting.set_TMD_size_exect==1
                    topologyandseqlog(count).(predictors{p})=tempseq(idx(1):(idx(1)+setting.set_TMD_size-1));
                elseif setting.set_TMD_size>tmd_size&&setting.set_TMD_size_exect==0
                    topologyandseqlog(count).(predictors{p})=tempseq(idx(1):idx(2));
                end
            elseif numberOfRegions==0
            end
        catch
        end
       end
    end
    textwaitbar(count,height(topologyandseqq), 'segmenting topology data')
    
end

qq=ones(1,length(predictors));
for pp=1:length(predictors)
    for K=1:height(topologyandseqq)
        try
            if ~isempty(topologyandseqlog(K).(predictors{pp}))
                seqdid=char(topologyandseqlog(K).(predictors{pp}));
                out.(predictors{pp})(qq(pp)).name=topologyandseqlog(K).name;
                out.(predictors{pp})(qq(pp)).hitname=topologyandseqlog(K).hitname;
                out.(predictors{pp})(qq(pp)).hitratio=topologyandseqlog(K).hitratio;
                out.(predictors{pp})(qq(pp)).tmdseq=seqdid;
                 out.(predictors{pp})(qq(pp)).TMDnum=topologyandseqlog(K).TMDnum;
                for q=1:length(setting.func)
                    if setting.get_functions_mean==1
                        out.(predictors{pp})(qq(pp)).(char(setting.func{q}))=mean(setting.func{q}(seqdid));
                    elseif setting.get_functions_mean==0
                        
                        if ~(setting.spesific_AA==0)
                            out.(predictors{pp})(qq(pp)).(char(setting.func{q}))=(setting.func{q}(seqdid(setting.spesific_AA)));
                        else
                            out.(predictors{pp})(qq(pp)).(char(setting.func{q}))=(setting.func{q}(seqdid));
                        end
                    end
                    
                end
                qq(pp)=qq(pp)+1;
            end
        catch
        end
        
    end
    textwaitbar(pp,length(predictors), 'calcolating TMD propertis')
    
end
if setting.get_functions_mean==1&&setting.save_sequance_analesis_data==1
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data'));

        for pp=1:length(predictors)
                tablename=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,'mean_analesis','_predection_',char(predictors{pp}),'_out.csv');
                writetable(struct2table(out.(predictors{pp})),tablename);
         end
else
        
used_predictors=fieldnames(out);
lis=ones(1,length(used_predictors));

for q2=1:length(setting.func)
    for pp2=1:length(used_predictors)
        
        outsplit.(group1).(used_predictors{pp2}).(char(setting.func{q2}))=[out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]==group1).(char(setting.func{q2}))]';
        outsplit.(group1).(used_predictors{pp2}).seq={out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]==group1).tmdseq};
        outsplit.(group2).(used_predictors{pp2}).(char(setting.func{q2}))=[out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]==group2).(char(setting.func{q2}))]';
        outsplit.(group2).(used_predictors{pp2}).seq={out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]==group2).tmdseq};
        outsplit.same.(used_predictors{pp2}).seq={out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]=='same').tmdseq};
        outsplit.same.(used_predictors{pp2}).(char(setting.func{q2}))=[out.(used_predictors{pp2})([out.(used_predictors{pp2}).hitname]=='same').(char(setting.func{q2}))]';
        lis(pp2)=lis(pp2)+1;
        
    end
end

lisa=ones(1,length(used_predictors));
textwaitbar(0,length(setting.func), 'looking for unique amino acids locos')
for q2=1:length(setting.func)
    for pp2=1:length(used_predictors)
        for ll=1:min(cellfun(@length,{out.(used_predictors{pp2}).tmdseq}))
            xi=outsplit.(group1).(used_predictors{pp2}).(char(setting.func{q2}))(:,ll);
            yi=outsplit.(group2).(used_predictors{pp2}).(char(setting.func{q2}))(:,ll);
            zi=outsplit.same.(used_predictors{pp2}).(char(setting.func{q2}))(:,ll);
            
            [h1,~] = ttest2(xi,yi,'Vartype','unequal');
            [h2,~] = ttest2(xi,zi,'Vartype','unequal');
            [h3,~] = ttest2(zi,yi,'Vartype','unequal');
            
            if h1&&h2&&~h3||h1&&~h2&&h3||h1&&h2&&h3||h1
                aalist.(used_predictors{pp2})(lisa(pp2)).hit=ll;
                aalist.(used_predictors{pp2})(lisa(pp2)).predictor=char(setting.func{q2});
                if h1&&h2&&~h3
                    aalist.(used_predictors{pp2})(lisa(pp2)).sig=group1;
                elseif h1&&~h2&&h3
                    aalist.(used_predictors{pp2})(lisa(pp2)).sig=group2;
                else
                    aalist.(used_predictors{pp2})(lisa(pp2)).sig='non group spesific';
                end
                lisa(pp2)=lisa(pp2)+1;
                
            end
            
            try
                %if setting.hits_uni==1
                    intrestingaa.(used_predictors{pp2}) = unique([unique([aalist.(used_predictors{pp2})(cellfun(@string,{aalist.(used_predictors{pp2}).sig}')==group1).hit]),unique([aalist.(used_predictors{pp2})(cellfun(@string,{aalist.(used_predictors{pp2}).sig}')==group2).hit])]);
                %elseif setting.hits_uni==0
                %    intrestingaa.(used_predictors{pp2}) = intersect(unique([aalist.(used_predictors{pp2})(cellfun(@string,{aalist.(used_predictors{pp2}).sig}')==group1).hit]),unique([aalist.(used_predictors{pp2})(cellfun(@string,{aalist.(used_predictors{pp2}).sig}')==group2).hit]));
                %end
            catch
            end
        end
    end
    textwaitbar(q2,length(setting.func), 'looking for unique amino acids locos')
end

    if setting.save_sequance_analesis_data==1
        
        mkdir(strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data'));
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
        textwaitbar(0,length(used_predictors), 'ploting')
        for pp=1:length(used_predictors)
                tablename=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_predection_',char(used_predictors{pp}),'_out.csv');
                writetable(struct2table(out.(used_predictors{pp})),tablename);
                clear hitlist hitlistsize
                %hitlist{3}=unique([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')==group1).hit]);                
                %hitlist{2}=unique([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')==group2).hit]);
                %hitlist{1}=unique([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')=='non group spesific').hit]);
                try
               % for lll=1:size(hitlist{1},2)
               % hitlistsize{1}(lll)=sum([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')=='non group spesific').hit]==hitlist{1}(lll));
               % end
                catch
                hitlistsize{1}=[];
                end
                try
               % for lll=1:size(hitlist{2},2)
               % hitlistsize{2}(lll)=sum([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')==group2).hit]==hitlist{2}(lll));
               % end
                catch
                hitlistsize{2}=[];
                end
                try
                for lll=1:size(hitlist{3},2)
                hitlistsize{3}(lll)=sum([aalist.(used_predictors{pp})(cellfun(@string,{aalist.(used_predictors{pp}).sig}')==group1).hit]==hitlist{3}(lll));
                end
                catch
                hitlistsize{3}=[];
                end
                if setting.do_seqlogo==1||setting.do_hlixwhell==1
                    
                FILENAME1=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_for_',group1,'_predection_',char(used_predictors{pp}),'_seq_logo.tif');
                FILENAME2=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_for_',group2,'_predection_',char(used_predictors{pp}),'_seq_logo.tif');
                FILENAME3=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_for_','same_predection_',char(used_predictors{pp}),'_seq_logo.tif');
                FILENAME_anotated=strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data',filesep,typeof,'_anotated_',char(used_predictors{pp}),'_seq_logo.tif');
                [w1, hfig1] =seqlogo_norm(outsplit.(group1).(used_predictors{pp}).seq,'Alphabet', 'AA');
                [w2, hfig2] =seqlogo_norm(outsplit.(group2).(used_predictors{pp}).seq,'Alphabet', 'AA');
                [w3, hfig3] =seqlogo_norm(outsplit.same.(used_predictors{pp}).seq,'Alphabet', 'AA');
                hfig1fraim=  frame2im(getframe(hfig1));
                hfig2fraim=  frame2im(getframe(hfig2));
                hfig3fraim=  frame2im(getframe(hfig3));
                 if setting.save_intermidiet_seqlogo==1
                     saveas(hfig1,FILENAME1)
                    saveas(hfig2,FILENAME2)
                    saveas(hfig3,FILENAME3)
                 end
                    close(hfig1)
                    close(hfig2)
                    close(hfig3)
                  
                if setting.do_seqlogo==1
                    plotinganotation(hfig1fraim,hfig2fraim,hfig3fraim,used_predictors{pp},group1,group2,hitlist,FILENAME_anotated,setting,hitlistsize)
                    close
                end



                 if setting.do_hlixwhell==1
                     
                 w{1,1}=w1{1,1};
                 w{1,2}=w1{1,2};
                 w{2,1}=w2{1,1};
                 w{2,2}=w2{1,2};
                 w{3,1}=w3{1,1};
                 w{3,2}=w3{1,2};
                 
                 for nn=1:3
                     missingaa=setdiff(setting.aa,w{nn,1});
                     if ~isempty(missingaa)
                         w{nn,1}=cat(1,w{nn,1},missingaa);
                         w{nn,2}=cat(1,w{nn,2},zeros(length(missingaa),size(w{nn,2},2)));
                     end
                 end
                 group={group1,group2,'same'};
                 helixwhell(w,used_predictors{pp},group,strcat(folder_path,filesep,'analesis_files',filesep,'sequance_analesis_data'),setting,folder_path);
                 end
                end
             textwaitbar(pp,length(used_predictors), 'ploting')
        end
    end
end
    disp('done!')
    
end