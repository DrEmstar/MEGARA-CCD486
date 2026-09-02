%MAKE_NEW_THARDAT Interactively rebuild the ThAr line-position calibration.
% Maintenance utility, not part of a normal reduction. Edit the input ThAr
% FITS filename and output filename below, run from a directory containing
% thar.dat, inspect the overlay, then select the corrections requested by
% configure_for_th. The output MAT file can be reviewed before installation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now lets try for a new thar image
    dataj = hercules_fitsread('J0889013.fit');
    %there is no blue_data_chop_value to be considered at this point
    blue_data_chop_value=0;
    figure(33)
    image(dataj,'cdatamapping','direct');colormap gray;%caxis([100 3000])
    hold on
    
    fid=fopen('thar.dat','r');
    C=textscan(fid,'%.0f %.4f %.4f %.4f %s %.3f %.3f %.3f %.3f %.3f %.3f');
    order=C{1};air=C{2};vac=C{3};wavnum=C{4};species=C{5};ul=C{6};uc=C{7};ur=C{8};vl=C{9};vc=C{10};vr=C{11};
    fclose(fid);
    % standard corrections to the text file
    uc=uc/0.015+36;vc=vc/0.015-200;
    ul=ul/0.015+36;vl=vl/0.015-200;
    ur=ur/0.015+36;vr=vr/0.015-200;
    %check start point looks similar to current image:
    plot(uc,vc,'r+')
    %shifts are much bigger in some parts but manageable, lets try the next step:
    [ufix,vfix] = configure_for_th(dataj,uc,vc,blue_data_chop_value);
    %the solution is much better! There are still some parts at the
    %extremes of the image that are not perfect - it might be an idea to
    %redo the thar distribution and list again, but not absolutely
    %necessary yet ...
    %now apply suggested fixes to uc and vc etc. 
    uc=uc+ufix;vc=vc+vfix;
    ul=ul+ufix;vl=vl+vfix;
    ur=ur+ufix;vr=vr+vfix;
    
    %and now save the result into a thardat.mat type file for use in reduction
save('thardat_p.mat','air','order','species','uc','ul','ur','vc','vl','vr','vac','wavnum')
    
