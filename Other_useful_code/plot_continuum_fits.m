%Plot the continuum fits from a completed MEGARA 2.1 reduction.

clearvars
clc

%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config.objectname='hd_39060';
config.plot_orders=[];       %Leave empty to plot every fitted order.
config.order_step=1;         %Use 5 or 10 for a quicker overview.
config.orders_per_figure=4;
config.maximum_spectra=50;   %Representative subset used for the median.

%These are normally read from the saved automatic continuum fit. Set any
%missing values here when reviewing a template fit with no automatic fit.
config.teff=[];
config.logg=[];
config.vsini=[];

%CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reduced_data_pointer=which('pointer_reduced_data_directory.m');
if isempty(reduced_data_pointer)
    error(['Could not find pointer_reduced_data_directory.m on the ' ...
        'MATLAB path.'])
end

reduced_data_directory=fileparts(reduced_data_pointer);
target_directory=fullfile(reduced_data_directory,config.objectname);
if ~isfolder(target_directory)
    error(['Could not find the reduced-data directory: ' target_directory])
end
reduced_frames_directory=fullfile(target_directory,'reduced frames');
supplementary_directory=fullfile(target_directory,'supplementary_data');

final_file=find_final_data_file(target_directory,config.objectname);
loaded_final=load(final_file,'finaldata');
if ~isfield(loaded_final,'finaldata')
    error(['The final-data file does not contain finaldata: ' final_file])
end
finaldata=loaded_final.finaldata;
continuum_metadata=finaldata;
diagnostics_file=fullfile(supplementary_directory, ...
    ['diagnostics_' config.objectname '.mat']);
if isfile(diagnostics_file)
    loaded_diagnostics=load(diagnostics_file,'diagnostics');
    if isfield(loaded_diagnostics,'diagnostics') && ...
            isfield(loaded_diagnostics.diagnostics,'continuum_fitting')
        continuum_metadata=loaded_diagnostics.diagnostics.continuum_fitting;
    end
end

[fit_file,fit_type]=find_continuum_fit_file( ...
    supplementary_directory,config.objectname,continuum_metadata);
loaded_fit=load(fit_file);
if ~isfield(loaded_fit,'datacf')
    error(['The continuum-fit file does not contain datacf: ' fit_file])
end
datacf=loaded_fit.datacf;

[teff,logg,vsini]=get_synthetic_parameters( ...
    loaded_fit,supplementary_directory,config);

order_numbers=find_fitted_orders(datacf);
if ~isempty(config.plot_orders)
    order_numbers=intersect(order_numbers,config.plot_orders,'stable');
end
order_numbers=order_numbers(1:config.order_step:end);
if isempty(order_numbers)
    error('No requested continuum-fit orders were found.')
end

reduced_files=dir(fullfile(reduced_frames_directory,'J*.mat'));
if isempty(reduced_files)
    error(['No reduced J*.mat spectra were found in ' reduced_frames_directory])
end
sample_indices=unique(round(linspace(1,numel(reduced_files), ...
    min(config.maximum_spectra,numel(reduced_files)))));
reduced_files=reduced_files(sample_indices);

fprintf('Using %d of %d reduced spectra to make the order medians\n', ...
    numel(reduced_files),numel(dir(fullfile(reduced_frames_directory,'J*.mat'))))
representatives=make_order_representatives( ...
    reduced_files,reduced_frames_directory,datacf,order_numbers);

wstart=floor(min(cellfun(@(x) min(x),representatives.wave)))-5;
wstop=ceil(max(cellfun(@(x) max(x),representatives.wave)))+5;
[synthetic_wave,synthetic_intensity]=make_synth_spectrum( ...
    teff,logg,vsini,0.5,wstart,wstop);

plot_order_fits(representatives,synthetic_wave,synthetic_intensity, ...
    config,fit_type,teff,logg,vsini)

fprintf('Displayed %d continuum-fit orders from %s\n', ...
    numel(order_numbers),fit_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function final_file=find_final_data_file(target_directory,objectname)

standard_file=fullfile(target_directory,['final_data_' objectname '.mat']);
extended_file=fullfile(target_directory, ...
    ['final_data_extended_' objectname '.mat']);
if isfile(standard_file)
    final_file=standard_file;
elseif isfile(extended_file)
    final_file=extended_file;
else
    error(['Could not find a final-data file for ' objectname '.'])
end

end

function [fit_file,fit_type]=find_continuum_fit_file( ...
    supplementary_directory,objectname,finaldata)

fit_directory=fullfile(supplementary_directory,'automatic_continuum_results');
if isfield(finaldata,'template_fit_method')
    fit_type='template';
    fit_file=fullfile(fit_directory, ...
        [objectname '_interpolated_spt_template_v1_standard.mat']);
elseif isfield(finaldata,'automatic_fit_method')
    fit_type='automatic';
    if isfield(finaldata,'extended_wavelength_range') && ...
            finaldata.extended_wavelength_range
        range_label='extended';
    else
        range_label='standard';
    end
    fit_file=fullfile(fit_directory, ...
        [objectname '_coarse_continuum_v6_' range_label '.mat']);
else
    fit_type='manual';
    fit_file=fullfile(supplementary_directory,['datacf_' objectname '.mat']);
    if ~isfile(fit_file)
        fit_file=fullfile(supplementary_directory,['datacf_' objectname]);
    end
end
if ~isfile(fit_file)
    error(['Could not find the saved ' fit_type ...
        ' continuum fit: ' fit_file])
end

end

function [teff,logg,vsini]=get_synthetic_parameters( ...
    loaded_fit,supplementary_directory,config)

teff=choose_parameter(config.teff,loaded_fit,'teff');
logg=choose_parameter(config.logg,loaded_fit,'logg');
vsini=choose_parameter(config.vsini,loaded_fit,'vsini');

if isempty(teff) || isempty(logg) || isempty(vsini)
    automatic_files=dir(fullfile(supplementary_directory, ...
        'automatic_continuum_results','*_coarse_continuum_v6_*.mat'));
    if ~isempty(automatic_files)
        [~,newest]=max([automatic_files.datenum]);
        previous_fit=load(fullfile(automatic_files(newest).folder, ...
            automatic_files(newest).name),'teff','logg','vsini');
        if isempty(teff), teff=choose_parameter([],previous_fit,'teff'); end
        if isempty(logg), logg=choose_parameter([],previous_fit,'logg'); end
        if isempty(vsini), vsini=choose_parameter([],previous_fit,'vsini'); end
    end
end

if isempty(teff) || isempty(logg)
    error(['The synthetic reference requires teff and logg. Set the ' ...
        'missing values in the USER INPUTS section.'])
end
if isempty(vsini) || ~isfinite(vsini)
    warning('No finite vsini was found. Using 0 km/s for the synthetic reference.')
    vsini=0;
end

end

function value=choose_parameter(user_value,loaded,field_name)

value=user_value;
if isempty(value) && isfield(loaded,field_name)
    candidate=loaded.(field_name);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value=candidate;
    end
end

end

function order_numbers=find_fitted_orders(datacf)

names=fieldnames(datacf);
order_numbers=[];
for field_index=1:numel(names)
    value=regexp(names{field_index},'^w(\d+)$','tokens','once');
    if ~isempty(value) && isfield(datacf,['int' value{1}])
        order_numbers(end+1)=str2double(value{1}); %#ok<AGROW>
    end
end
order_numbers=sort(order_numbers);

end

function representatives=make_order_representatives( ...
    reduced_files,reduced_frames_directory,datacf,order_numbers)

number_orders=numel(order_numbers);
number_spectra=numel(reduced_files);
representatives.order=order_numbers;
representatives.wave=cell(1,number_orders);
samples=cell(1,number_orders);
for order_index=1:number_orders
    wave=datacf.(['w' num2str(order_numbers(order_index))]);
    representatives.wave{order_index}=wave(:).';
    samples{order_index}=nan(number_spectra,numel(wave));
end

for spectrum_index=1:number_spectra
    loaded=load(fullfile(reduced_frames_directory,reduced_files(spectrum_index).name));
    reduced_frame=get_primary_struct(loaded,reduced_files(spectrum_index).name);
    for order_index=1:number_orders
        order_number=order_numbers(order_index);
        order_field=['order_' num2str(order_number)];
        if ~isfield(reduced_frame,order_field), continue; end
        order_data=reduced_frame.(order_field);
        if size(order_data,2)<2, continue; end
        valid=isfinite(order_data(:,1)) & isfinite(order_data(:,2));
        if sum(valid)<2, continue; end
        samples{order_index}(spectrum_index,:)=interp1( ...
            order_data(valid,1),order_data(valid,2), ...
            representatives.wave{order_index},'linear',NaN);
    end
end

representatives.intensity=cell(1,number_orders);
representatives.continuum=cell(1,number_orders);
for order_index=1:number_orders
    order_number=order_numbers(order_index);
    representatives.intensity{order_index}= ...
        median(samples{order_index},1,'omitnan');
    representatives.continuum{order_index}= ...
        datacf.(['int' num2str(order_number)])(:).';
end

end

function reduced_frame=get_primary_struct(loaded,file_name)

names=fieldnames(loaded);
reduced_frame=[];
for field_index=1:numel(names)
    if isstruct(loaded.(names{field_index}))
        reduced_frame=loaded.(names{field_index});
        break
    end
end
if isempty(reduced_frame)
    error(['Could not find the reduced-frame structure in ' file_name])
end

end

function plot_order_fits(representatives,synthetic_wave, ...
    synthetic_intensity,config,fit_type,teff,logg,vsini)

number_orders=numel(representatives.order);
number_figures=ceil(number_orders/config.orders_per_figure);
for figure_index=1:number_figures
    first_index=(figure_index-1)*config.orders_per_figure+1;
    last_index=min(first_index+config.orders_per_figure-1,number_orders);
    displayed_indices=first_index:last_index;

    figure('Color','w','Name',sprintf('%s continuum fits %d of %d', ...
        config.objectname,figure_index,number_figures))
    layout=tiledlayout(numel(displayed_indices),2, ...
        'TileSpacing','compact','Padding','compact');
    title(layout,sprintf(['%s: %s continuum, Teff %.0f K, ' ...
        'logg %.1f, vsini %.1f km/s'],strrep(config.objectname,'_',' '), ...
        fit_type,teff,logg,vsini))

    for order_index=displayed_indices
        order_number=representatives.order(order_index);
        wave=representatives.wave{order_index};
        observed=representatives.intensity{order_index};
        continuum=representatives.continuum{order_index};
        synthetic=interp1(synthetic_wave,synthetic_intensity,wave, ...
            'linear',NaN);
        valid_continuum=isfinite(continuum) & continuum>0;
        normalised=observed./continuum;
        normalised(~valid_continuum)=NaN;

        nexttile
        plot(wave,observed,'Color',[0.45 0.45 0.45],'LineWidth',0.5)
        hold on
        plot(wave,continuum,'r','LineWidth',1.2)
        ylabel('Intensity')
        title(sprintf('Order %d fit',order_number))
        legend('Observed median','Continuum fit','Location','best')
        grid on

        nexttile
        plot(wave,normalised,'k','LineWidth',0.6)
        hold on
        plot(wave,synthetic,'b','LineWidth',0.7)
        yline(1,':','Color',[0.5 0.5 0.5])
        ylabel('Normalised intensity')
        title(sprintf('Order %d normalised',order_number))
        legend('Observed / fit','Synthetic reference','Location','best')
        grid on
    end
    xlabel(layout,'Wavelength (Angstrom)')
end

end
