function [synthwave,synthint,wavelengths,ele,elem,widths]=make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop)

validateattributes(teff,{'numeric'},{'scalar','real','finite','positive'},mfilename,'teff')
validateattributes(logg,{'numeric'},{'scalar','real','finite'},mfilename,'logg')
validateattributes(vsini,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'vsini')
validateattributes(sigm,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'sigm')
validateattributes(wstart,{'numeric'},{'scalar','real','finite','positive'},mfilename,'wstart')
validateattributes(wstop,{'numeric'},{'scalar','real','finite','positive'},mfilename,'wstop')

if wstop<=wstart
    error('wstop must be greater than wstart.')
end

original_directory=pwd;
synspec_pointer_path=which('synspec_pointer.m');

if isempty(synspec_pointer_path)
    error('Could not find synspec_pointer.m on the MATLAB path.')
end

synspec_directory=fileparts(synspec_pointer_path);
cd(synspec_directory)
synspec_cleanup=onCleanup(@() finish_synspec_run(original_directory,synspec_directory)); %#ok<NASGU>

[model_teff,model_logg]=select_nearest_atmosphere_model( ...
    teff,logg,synspec_directory);

if model_teff~=teff || model_logg~=logg
    fprintf(['Using nearest available atmosphere model: Teff %.0f K, ' ...
        'logg %.1f\n'],model_teff,model_logg)
end

write_new_file_for_synspec(model_teff,model_logg);
write_mod55_file_for_synspec(wstart,wstop);

remove_stale_synspec_working_files()
delete_file_if_present('mod.7')

[rmod_status,rmod_message]=system('./rmod');

if rmod_status~=0
    error('rmod failed with status %d.\n%s',rmod_status,rmod_message)
end

check_nonempty_output_file('mod.7','rmod',rmod_message)

remove_stale_synspec_working_files()
delete_file_if_present('output.7')
delete_file_if_present('output.12')
delete_file_if_present('output.17')

[synspec_status,synspec_message]=system('./KSynspec output mod gf3000');

check_nonempty_output_file('output.7','KSynspec',synspec_message)
check_nonempty_output_file('output.12','KSynspec',synspec_message)
check_nonempty_output_file('output.17','KSynspec',synspec_message)

%KSynspec can return status 1 when optional output files are absent.

[wavelengths,ele,elem,widths]=read_in_output12;

if isempty(wavelengths) || isempty(ele) || isempty(elem) || isempty(widths)
    error('SynSpec output.12 was read but did not contain a usable line list.')
end

write_r_file_for_rotin3(vsini,sigm,wstart,wstop);
delete_file_if_present('output.11')

[rotin3_status,rotin3_message]=system('./rotin3 < r.dat');

if rotin3_status~=0
    error('rotin3 failed with status %d.\n%s',rotin3_status,rotin3_message)
end

check_nonempty_output_file('output.11','rotin3',rotin3_message)

[synthwave,synthint]=read_in_output11;

if isempty(synthwave) || isempty(synthint)
    error('rotin3 output.11 was read but did not contain a usable synthetic spectrum.')
end

if numel(synthwave)~=numel(synthint)
    error('Synthetic wavelength and intensity vectors have different lengths.')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model_teff,model_logg]=select_nearest_atmosphere_model( ...
    teff,logg,synspec_directory)

model_file=fullfile(synspec_directory,'models','ap00k2odfnew.dat');
if ~isfile(model_file)
    error(['Could not find the SynSpec atmosphere-model grid: ' model_file])
end

file_identifier=fopen(model_file,'r');
if file_identifier==-1
    error(['Could not open the SynSpec atmosphere-model grid: ' model_file])
end
cleanup_file=onCleanup(@() fclose(file_identifier)); %#ok<NASGU>

available_teff=[];
available_logg=[];
line=fgetl(file_identifier);
while ischar(line)
    values=regexp(line, ...
        '^\s*TEFF\s+([0-9.]+)\s+GRAVITY\s+([-+0-9.]+)', ...
        'tokens','once');
    if ~isempty(values)
        available_teff(end+1)=str2double(values{1}); %#ok<AGROW>
        available_logg(end+1)=str2double(values{2}); %#ok<AGROW>
    end
    line=fgetl(file_identifier);
end

valid=isfinite(available_teff) & isfinite(available_logg);
available_teff=available_teff(valid);
available_logg=available_logg(valid);
if isempty(available_teff)
    error(['No atmosphere-model headers were found in ' model_file])
end

temperature_values=unique(available_teff);
[~,temperature_index]=min(abs(temperature_values-teff));
model_teff=temperature_values(temperature_index);

at_temperature=available_teff==model_teff;
gravity_values=unique(available_logg(at_temperature));
[~,gravity_index]=min(abs(gravity_values-logg));
model_logg=gravity_values(gravity_index);

end

function finish_synspec_run(original_directory,synspec_directory)

try
    if isfolder(synspec_directory)
        cd(synspec_directory)
        remove_stale_synspec_working_files()
    end
catch cleanup_error
    warning('Could not completely clean the SynSpec working directory: %s',cleanup_error.message)
end

cd(original_directory)

end

function remove_stale_synspec_working_files()

working_files=[dir('fort.*');dir('f99')];

for file_index=1:numel(working_files)

    if working_files(file_index).isdir
        continue
    end

    file_path=fullfile(working_files(file_index).folder,working_files(file_index).name);

    try
        delete(file_path)
    catch cleanup_error
        error('Could not remove stale SynSpec working file %s: %s',file_path,cleanup_error.message)
    end

end

end

function delete_file_if_present(file_name)

file_information=dir(file_name);

for file_index=1:numel(file_information)

    if file_information(file_index).isdir
        continue
    end

    file_path=fullfile(file_information(file_index).folder,file_information(file_index).name);
    delete(file_path)

end

end

function check_nonempty_output_file(file_name,program_name,program_message)

file_information=dir(file_name);

if numel(file_information)~=1 || file_information.isdir
    error('%s did not create %s.\n%s',program_name,file_name,program_message)
end

if file_information.bytes<=0
    error('%s created an empty %s file.\n%s',program_name,file_name,program_message)
end

end
