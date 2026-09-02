function datacf=prepare_spt_template_continuum_fit(data2,objectname,teff, ...
    firstorder,lastorder,extended_wavelength_range)
%Interpolate the SpT continuum grid and apply it directly to each order.

if isempty(teff) || ~isnumeric(teff) || ~isscalar(teff) || ~isfinite(teff)
    error('SpT template continuum fitting requires a finite scalar teff.')
end
if extended_wavelength_range
    error(['The SpT continuum grid supports only the standard 3869-6866 ' ...
        'Angstrom wavelength range.'])
end

grid_file=find_spt_template_grid_file();
loaded_grid=load(grid_file,'template_grid');
if ~isfield(loaded_grid,'template_grid')
    error(['Template-grid file does not contain template_grid: ' grid_file])
end
template_grid=loaded_grid.template_grid;
required_fields={'teff','teff_lower','teff_upper','star'};
for field_index=1:numel(required_fields)
    if ~isfield(template_grid,required_fields{field_index})
        error(['template_grid does not contain ' required_fields{field_index} '.'])
    end
end

grid_teff=template_grid.teff(:);
if any(diff(grid_teff)<=0)
    error('Template-grid temperatures must be strictly increasing.')
end
minimum_supported_temperature=min(template_grid.teff_lower);
maximum_supported_temperature=max(template_grid.teff_upper);
if teff<minimum_supported_temperature || teff>maximum_supported_temperature
    error(['teff ' num2str(teff) ' K is outside the SpT template range ' ...
        num2str(minimum_supported_temperature) '-' ...
        num2str(maximum_supported_temperature) ' K.'])
end

[lower_index,upper_index,temperature_weight]= ...
    find_temperature_bracket(grid_teff,teff);
lower_star=template_grid.star{lower_index};
upper_star=template_grid.star{upper_index};
if lower_index==upper_index
    fprintf('Using SpT continuum template %s at %.0f K\n',lower_star, ...
        grid_teff(lower_index))
else
    fprintf(['Interpolating SpT continuum templates %s (%.0f K) and ' ...
        '%s (%.0f K)\n'],lower_star,grid_teff(lower_index), ...
        upper_star,grid_teff(upper_index))
end

config.minimum_template_coverage=0.80;
config.maximum_template_shape_range=0.50;
config.minimum_scale_points=30;
config.scale_percentile=99.5;

datacf=struct();
fit_summary=cell(0,8);
for order_number=firstorder:lastorder
    wave_field=['wave' num2str(order_number)];
    intensity_field=['int' num2str(order_number)];
    template_wave_field=['w' num2str(order_number)];
    if ~isfield(data2,wave_field) || ~isfield(data2,intensity_field)
        continue
    end

    wave=data2.(wave_field)(:).';
    continuum=nan(size(wave));
    coverage=0;
    shape_range=NaN;
    scale=NaN;
    status='unavailable';

    if isfield(template_grid,template_wave_field) && ...
            isfield(template_grid,intensity_field)
        grid_wave=template_grid.(template_wave_field)(:).';
        continuum_grid=template_grid.(intensity_field);
        lower_shape=interpolate_template_row(grid_wave, ...
            continuum_grid(lower_index,:),wave);
        if upper_index==lower_index
            template_shape=lower_shape;
        else
            upper_shape=interpolate_template_row(grid_wave, ...
                continuum_grid(upper_index,:),wave);
            template_shape=(1-temperature_weight).*lower_shape+ ...
                temperature_weight.*upper_shape;
        end

        valid_shape=isfinite(template_shape) & template_shape>0;
        coverage=mean(valid_shape);
        if any(valid_shape)
            shape_limits=prctile(template_shape(valid_shape),[5 95]);
            shape_range=diff(shape_limits);
        end

        representative=median(data2.(intensity_field),1,'omitnan');
        valid_scale=valid_shape & isfinite(representative) & representative>0;
        if coverage>=config.minimum_template_coverage && ...
                isfinite(shape_range) && ...
                shape_range<=config.maximum_template_shape_range && ...
                sum(valid_scale)>=config.minimum_scale_points
            scale=prctile(representative(valid_scale)./ ...
                template_shape(valid_scale),config.scale_percentile);
            if isfinite(scale) && scale>0
                continuum=template_shape.*scale;
                status='interpolated';
            end
        end
    end

    if strcmp(status,'unavailable')
        fprintf(['Template order %.0f is unavailable or unstable at %.0f K; ' ...
            'excluding it from the merge\n'],order_number,teff)
    end
    datacf.(template_wave_field)=wave;
    datacf.(intensity_field)=continuum;
    fit_summary(end+1,:)={order_number,lower_star,upper_star, ...
        temperature_weight,coverage,shape_range,scale,status}; %#ok<AGROW>
end

datacf.template_fit_method='spt_grid_linear_temperature_v1';
datacf.template_fit_summary=cell2table(fit_summary,'VariableNames', ...
    {'Order','LowerStar','UpperStar','TemperatureWeight','Coverage', ...
    'ShapeRange','Scale','Status'});
datacf.template_grid_file=grid_file;
datacf.template_teff=teff;

output_directory=fullfile(pwd,'supplementary_data', ...
    'automatic_continuum_results');
if ~isfolder(output_directory), mkdir(output_directory); end
output_file=fullfile(output_directory, ...
    [objectname '_interpolated_spt_template_v1_standard.mat']);
save(output_file,'datacf','config','teff','grid_file','-v7.3')

end

function grid_file=find_spt_template_grid_file()

grid_name='spt_continuum_template_grid_v1_standard.mat';
grid_file=which(grid_name);
if ~isempty(grid_file), return; end

post_reduction_directory=fileparts(mfilename('fullpath'));
megara_directory=fileparts(fileparts(post_reduction_directory));
grid_file=fullfile(megara_directory,'Other_useful_code',grid_name);
if ~isfile(grid_file)
    error(['Could not find ' grid_name '. Build the SpT template grid ' ...
        'before using options.template_continuum_fit.'])
end

end

function [lower_index,upper_index,weight]= ...
    find_temperature_bracket(grid_teff,teff)

[nearest_offset,nearest_index]=min(abs(grid_teff-teff));
if nearest_offset<1e-8
    lower_index=nearest_index;
    upper_index=nearest_index;
    weight=0;
elseif teff<grid_teff(1)
    lower_index=1;
    upper_index=1;
    weight=0;
elseif teff>grid_teff(end)
    lower_index=numel(grid_teff);
    upper_index=lower_index;
    weight=0;
else
    lower_index=find(grid_teff<teff,1,'last');
    upper_index=find(grid_teff>teff,1,'first');
    weight=(teff-grid_teff(lower_index))/ ...
        (grid_teff(upper_index)-grid_teff(lower_index));
end

end

function interpolated_shape=interpolate_template_row( ...
    grid_wave,template_shape,observed_wave)

interpolated_shape=nan(size(observed_wave));
valid=isfinite(grid_wave) & isfinite(template_shape) & template_shape>0;
if sum(valid)<2, return; end

valid_wave=grid_wave(valid);
valid_shape=template_shape(valid);
[valid_wave,unique_index]=unique(valid_wave,'stable');
valid_shape=valid_shape(unique_index);
if numel(valid_wave)<2, return; end
inside=observed_wave>=valid_wave(1) & observed_wave<=valid_wave(end);
interpolated_shape(inside)=interp1(valid_wave,valid_shape, ...
    observed_wave(inside),'linear');

end
