function signal_to_noise=get_signal_to_noise_of_image(filename,allorders,blue_data_chop,header2)
% this calculates the signal to noise at an approximate wavelength of 5777 Angstroms
try
    image=raw_hercules_fitsread(filename);
    if isempty(header2)
        header = get_header(filename);
        header2=get_required_headers_from_header(header);
    else
    end
    exp_time=header2.EXPTIME;
    extracted_stellar= extract_all_orders_no_background(allorders,image);
    stellar_data=redo_order_numbering_simple(extracted_stellar,allorders);
    stellar_order=stellar_data.order_100;
    signal_to_noise=compute_signal_to_noise_spectrum(stellar_order,exp_time);
catch
    signal_to_noise=0;
end