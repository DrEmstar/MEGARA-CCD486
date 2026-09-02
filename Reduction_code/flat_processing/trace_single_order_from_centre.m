function [pos,centres,widths]=trace_single_order_from_centre(summed_flat,starting_centre,starting_width)

%Validate inputs
validate_trace_single_order_inputs(summed_flat,starting_centre,starting_width);

%Set tracing parameters
step=10;
failtimes=20;

%Get image size
image_size=size(summed_flat);

%trace from centre to left in 'step' pixel steps
image_left_right_centre=round(image_size(2)/2); %centre pixel
left_positions=image_left_right_centre:-step:1;
if isempty(left_positions) || left_positions(end)~=1
    left_positions(end+1)=1;
end

%Trace left side
[left_positions,left_centres,left_widths]=trace_order_direction(summed_flat,left_positions,starting_centre,starting_width,'left',failtimes);

%now go centre to right
right_start=image_left_right_centre+step;
if right_start>image_size(2)
    right_positions=[];
else
    right_positions=right_start:step:image_size(2);
    if right_positions(end)~=image_size(2)
        right_positions(end+1)=image_size(2);
    end
end

%Trace right side
[right_positions,right_centres,right_widths]=trace_order_direction(summed_flat,right_positions,starting_centre,starting_width,'right',failtimes);

%arrange the fits
left_centres=flipud(left_centres);
left_widths=flipud(left_widths);
left_positions=flipud(left_positions);
pos=cat(1,left_positions,right_positions);
centres=cat(1,left_centres,right_centres);
widths=cat(1,left_widths,right_widths);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_trace_single_order_inputs(summed_flat,starting_centre,starting_width)

%Check image input
if ~isnumeric(summed_flat) || ndims(summed_flat)~=2, error('summed_flat must be a numeric 2D array.'); end
if any(~isfinite(summed_flat(:))), warning('summed_flat contains non-finite values.'); end

%Check starting centre
if ~isnumeric(starting_centre) || ~isscalar(starting_centre) || ~isfinite(starting_centre), error('starting_centre must be a finite scalar numeric value.'); end

%Check starting width
if ~isnumeric(starting_width) || ~isscalar(starting_width) || ~isfinite(starting_width) || starting_width<=0, error('starting_width must be a finite positive scalar numeric value.'); end

end

function [valid_positions,valid_centres,valid_widths]=trace_order_direction(summed_flat,trace_positions,starting_centre,starting_width,direction_name,failtimes)

%Initialise tracing arrays
image_size=size(summed_flat);
centres=zeros(numel(trace_positions),1);
widths=zeros(numel(trace_positions),1);
valid_trace=false(numel(trace_positions),1);
current_centre=starting_centre;
consecutive_failures=0;

%Trace order through requested columns
for position_index=1:numel(trace_positions)
    column_position=trace_positions(position_index);
    row_start=round(current_centre-3*starting_width);
    row_stop=round(current_centre+3*starting_width);
    if row_start<1 || row_stop>image_size(1) %then outside image so exit
        if strcmp(direction_name,'left')
            break
        else
            consecutive_failures=consecutive_failures+1;
            if consecutive_failures>=failtimes, break; end
            continue
        end
    end
    if strcmp(direction_name,'left')
        column_start=max(1,column_position);
        column_stop=min(image_size(2),column_position+3);
    else
        column_start=max(1,column_position-3);
        column_stop=min(image_size(2),column_position);
    end
    xdata=(row_start:row_stop)';
    ydata=mean(summed_flat(row_start:row_stop,column_start:column_stop),2,'omitnan');
    ydata=ydata(:)-min(ydata);
    [measured_centre,measured_width]=find_max_and_fwhm(xdata,ydata);
    if abs(measured_centre-current_centre)>5 || measured_width<2.5 || measured_width>11 || ~isfinite(measured_centre) || ~isfinite(measured_width)
        %if fail failtimes consecutive times then exit
        consecutive_failures=consecutive_failures+1;
        if consecutive_failures>=failtimes
            break
        end
        continue
    end
    consecutive_failures=0;
    centres(position_index)=measured_centre;
    widths(position_index)=measured_width;
    valid_trace(position_index)=true;
    current_centre=measured_centre;
end

%Return valid trace points
valid_positions=trace_positions(:);
valid_positions=valid_positions(valid_trace);
valid_centres=centres(valid_trace);
valid_widths=widths(valid_trace);

end