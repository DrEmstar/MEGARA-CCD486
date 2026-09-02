function display_directory_file_information()
list=dir('J*.fit');
disp(' ')
disp('------------------------------------------------------------')
disp('Filename      HERCEXPT  JDmid       FIBRE   Object S/N EXPTIME')
disp('------------------------------------------------------------')
fid=fopen('file_information.txt','w');
%Test if reduction has run
if isempty(dir('flat_info*'))
    no_sn=1;
else
    no_sn=0;
    try
        blue_data_chop=500;
        load flat_info_blue_500.mat
    catch
        no_sn=1;
    end
end
for s=1:numel(list)
    nam=list(s).name;nam(end-3:end)=[];
    eval(['header = get_header(''' list(s).name ''');'])
    temp = get_required_headers_from_header(header);
    temp=check_hercules_headers(temp);
    eval(['info.list{s}=''' nam ''';'])
    eval(['info.' nam '=temp;'])
    type=temp.HERCEXPT;
    
    if strcmp(type,'Stellar')
        if no_sn==1
            signal_to_noise=NaN;
        else
            signal_to_noise=get_signal_to_noise_of_image(list(s).name,allorders,blue_data_chop,temp);
        end
        fprintf('%s  %7s   %.4f    %d     %s     %2.1f     %7.2f\n',list(s).name,temp.HERCEXPT,temp.MJD_OBS,temp.HERCFIB,temp.OBJECT,signal_to_noise,temp.EXPTIME);
    else
        signal_to_noise='     ';
        fprintf('%s  %7s   %.4f    %d     %s     %s     %7.2f\n',list(s).name,temp.HERCEXPT,temp.MJD_OBS,temp.HERCFIB,temp.OBJECT,signal_to_noise,temp.EXPTIME);
    end
   
    fprintf(fid,'%s  %7s   %.4f    %d      %s     %7.2f\n',list(s).name,temp.HERCEXPT,temp.MJD_OBS,temp.HERCFIB,temp.OBJECT,temp.EXPTIME);
end
disp('---------------------------------------------')
fclose(fid);

