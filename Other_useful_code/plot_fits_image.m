function []=plot_fits_image(filename)
%This plots the frame for inspection. This code accepts inputs as filenames in .fits formats or pre-loaded data matrices in matlab.
if isa(filename,'double')
    image_data=filename;
    image(image_data,'cdatamapping','scaled');colormap gray;caxis([100 3000])
else
    image_data = raw_hercules_fitsread(filename);
    image(image_data,'cdatamapping','scaled');colormap gray;caxis([100 3000])
end





