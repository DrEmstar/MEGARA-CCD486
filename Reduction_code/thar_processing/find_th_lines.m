function [thfitinfo,numbad_thar_flag,move_file_flag]=find_th_lines(thdata,order,ul,uc,ur,plot_q)

%Initialise outputs
thfitinfo=[];
numbad_thar_flag=0;
move_file_flag=0;

%Validate plotting option
if isempty(plot_q) || numel(plot_q)~=1
    plot_q=0;
end

%Check input line-list lengths
number_of_lines=numel(uc);
if numel(order)~=number_of_lines || numel(ul)~=number_of_lines || numel(ur)~=number_of_lines, error('order, ul, uc, and ur must have the same number of elements.'); end

%Initialise counters
number_candidate_lines=0;
number_accepted_lines=0;

%Fit each ThAr line
for line_index=1:number_of_lines
    order_field_name=['order_' num2str(order(line_index))];
    if ~isfield(thdata,order_field_name), continue; end
    if ~isfield(thdata.(order_field_name),'xax') || ~isfield(thdata.(order_field_name),'summed_data'), continue; end
    xdata=thdata.(order_field_name).xax;
    ydata=thdata.(order_field_name).summed_data;
    xdata=xdata(:);
    ydata=ydata(:);
    if isempty(xdata) || isempty(ydata), continue; end
    if numel(xdata)~=numel(ydata), continue; end
    if ul(line_index)>=min(xdata) && ur(line_index)<=max(xdata)
        number_candidate_lines=number_candidate_lines+1;
        selected_pixels=xdata>=round(ul(line_index)) & xdata<=round(ur(line_index));
        xdata=xdata(selected_pixels);
        ydata=ydata(selected_pixels);
    else
        continue
    end
    valid_line_points=isfinite(xdata) & isfinite(ydata);
    xdata=xdata(valid_line_points);
    ydata=ydata(valid_line_points);
    if numel(xdata)<5, continue; end
    ydata=ydata-min(ydata);
    mean_ydata=mean(ydata,'omitnan');
    if isfinite(mean_ydata) && mean_ydata~=0
        ydata=ydata/mean_ydata;
    else
        continue
    end
    valid_line_points=isfinite(xdata) & isfinite(ydata);
    xdata=xdata(valid_line_points);
    ydata=ydata(valid_line_points);
    if numel(xdata)<5, continue; end
    lower_bounds=[0 ul(line_index) 1 -0.5 -1];
    startpoint=[1 uc(line_index) 4 0 0];
    upper_bounds=[50 ur(line_index) 15 1.0 1];
    [fit_parameters,resnorm,residual,exitflag,fitted]=fit_th_line(xdata,ydata,lower_bounds,upper_bounds,startpoint);
    if plot_q==1
        fprintf('resnorm=%.3f\n',resnorm)
        fprintf('hght=%.2f posn=%.1f FWHM=%.2f pos_ud=%.4f slope=%.4f\n',fit_parameters(1),fit_parameters(2),fit_parameters(3),fit_parameters(4),fit_parameters(5))
        disp(' ')
        plot(xdata,ydata)
        title(['line number ' num2str(line_index)])
        hold on
        plot(xdata,fitted,'r')
        hold off
        user_response=input('enter x to stop displaying ThAr fits -> ','s');
        if strcmp(user_response,'x')
            plot_q=0;
        end
    end
    if resnorm<=2 && fit_parameters(1)>0.8 && fit_parameters(3)>=1.5 && fit_parameters(3)<=9
        number_accepted_lines=number_accepted_lines+1;
        thfitinfo(number_accepted_lines,:)=cat(2,line_index,resnorm,fit_parameters);
    end
end

%Flag weak ThAr solutions
number_fitted_lines=size(thfitinfo,1);
if number_fitted_lines<800
    numbad_thar_flag=1;
else
    numbad_thar_flag=0;
end

%Report fitted lines
fprintf('%.0f lines found\n',number_fitted_lines)

%Reject unusable ThAr image
if number_fitted_lines<200
    fprintf('bad image definition- skipping\n\n')
    move_file_flag=1;
    return
end

end