function [result,thecfit]=continuum_fitting_inline(wave,int,synthint)
%This program allows manual editing of continuum anchor points.

%Initialise fit
cfit=[];
[thecfit,result]=update_plot(wave,int,cfit,synthint);

%Interactive edit loop
while true
    disp('1 - add point')
    disp('2 - delete point')
    disp('3 - move point')
    disp('4 - exit')
    disp(' ')
    entry=input('-> ');
    if isempty(entry)
        entry=0;
    end
    switch entry
        case 1
            [result,thecfit,cfit]=add_point(wave,int,cfit,synthint);
        case 2
            [result,thecfit,cfit]=delete_point(wave,int,cfit,synthint);
        case 3
            [result,thecfit,cfit]=move_point(wave,int,cfit,synthint);
        case 4
            break
        otherwise
            disp('not one of the available options ...')
    end
    disp(' ')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=add_point(wave,int,cfit,synthint)

%Add selected point
zoom off
pan off
figure(100)
disp('waiting for mouse button press within axes')
selected_point=ginput(1);
if isempty(selected_point)
    disp('no point selected')
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
    return
end
cfit=sortrows(cat(1,cfit,selected_point(:,1:2)),1);
cfit=remove_duplicate_wavelength_points(cfit);
try
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
catch
    disp('error in update plot, possibly two points at same place, try again!')
    thecfit=[];
    result=[];
    return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=delete_point(wave,int,cfit,synthint)

%Delete nearest point
if isempty(cfit)
    disp('no points to delete!')
    thecfit=[];
    result=[];
    return
end
zoom off
pan off
figure(100)
disp('waiting for mouse button press within axes')
xy_values=ginput(1);
if isempty(xy_values)
    disp('no point selected')
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
    return
end
xy_nearest=findnearest2D(cfit,xy_values(:,1:2));
cfit(xy_nearest,:)=[];
try
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
catch
    disp('error in update plot, possibly two points at same place, try again!')
    thecfit=[];
    result=[];
    return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=move_point(wave,int,cfit,synthint)

%Move nearest point
if isempty(cfit)
    disp('no points to move!')
    thecfit=[];
    result=[];
    return
end
zoom off
pan off
disp('waiting for mouse button press within axes')
figure(100)
while true
    button_press=waitforbuttonpress;
    if button_press==0
        pos=get(gca,'CurrentPoint');
        pos=pos(1,1:2);
        xy_nearest=findnearest2D(cfit,pos);
        plot(cfit(xy_nearest,1),cfit(xy_nearest,2),'ko')
        disp('circled point will move to position of next mouse button press ...')
        break
    end
end
disp('waiting for mouse button press within axes')
figure(100)
while true
    button_press=waitforbuttonpress;
    if button_press==0
        pos=get(gca,'CurrentPoint');
        pos=pos(1,1:2);
        cfit(xy_nearest,:)=pos;
        cfit=sortrows(cfit,1);
        cfit=remove_duplicate_wavelength_points(cfit);
        try
            [thecfit,result]=update_plot(wave,int,cfit,synthint);
        catch
            disp('error in update plot, possibly two points at same place, try again!')
            thecfit=[];
            result=[];
            return
        end
        break
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thecfit,result]=update_plot(wave,int,cfit,synthint)

%Update diagnostic plot
figure(100)
plot_axes(1)=subplot(2,1,1,'replace');
hold on
title('Input and continuum fit')
ylabel('Intensity (unknown)')
xlabel('Wavelength')
plot(wave,int,'Parent',plot_axes(1))
safe_synthint=synthint;
safe_synthint(safe_synthint==0)=NaN;
plot(wave,int./safe_synthint,'g','Parent',plot_axes(1))

%Plot selected continuum points
if size(cfit,1)>0
    plot(cfit(:,1),cfit(:,2),'+r','Parent',plot_axes(1))
end

%Fit continuum if enough points exist
if size(cfit,1)>=7
    thecfit=spline(cfit(:,1),cfit(:,2),wave);
    result=int./thecfit;
    plot(wave,thecfit,'g','Parent',plot_axes(1))
    plot_axes(2)=subplot(2,1,2,'replace');
    hold on
    title('Resulting fitted data')
    xlabel('Wavelength (Angstroms)')
    ylabel('Normalised Intensity')
    axis([-inf inf -inf inf])
    linkaxes(plot_axes,'x');
    plot(wave,synthint,'r','Parent',plot_axes(2))
    plot(wave,result,'Parent',plot_axes(2))
else
    thecfit=[];
    result=[];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy_nearest=findnearest2D(vector2D,xy_values)

%Find nearest 2D point
if isempty(vector2D)
    xy_nearest=[];
    return
end
sz=size(vector2D);
if sz(2)>2
    error('problem with input vector for findnearest2D')
end
dist=sqrt((vector2D(:,1)-xy_values(1)).^2+(vector2D(:,2)-xy_values(2)).^2);
xy_nearest=find(dist==min(dist));
if numel(xy_nearest)>1
    xy_nearest=xy_nearest(ceil(numel(xy_nearest)/2));
end

end

function cfit=remove_duplicate_wavelength_points(cfit)

%Remove duplicate x values
if isempty(cfit)
    return
end
[~,unique_indices]=unique(cfit(:,1),'stable');
if numel(unique_indices)<size(cfit,1)
    disp('duplicate wavelength point removed')
end
cfit=cfit(sort(unique_indices),:);

end