function bcorr=get_barycentric_correction(starname,expdate,midtime,filename)
%GET_BARYCENTRIC_CORRECTION Calculate the HERCULES correction in MATLAB.
% The text return value preserves the interface used by existing MEGARA
% callers. Catalogue targets include parallax, proper motion and RV. An
% unidentified target falls back to its FITS coordinates with a warning.

numeric_bcorr=calculate_barycentric_correction( ...
    expdate,midtime,starname);

if ~isfinite(numeric_bcorr)
    warning(['Target %s was not identified in the stellar catalogue. ' ...
        'Using FITS coordinates without parallax, proper motion or ' ...
        'radial velocity.'],char(string(starname)))
    [right_ascension,declination]=read_fallback_coordinates(filename);
    numeric_bcorr=calculate_barycentric_correction( ...
        expdate,midtime,right_ascension,declination);
end

if isfinite(numeric_bcorr)
    bcorr=sprintf('%.12f',numeric_bcorr);
else
    bcorr='NaN';
    warning('Barycentric correction was not calculated correctly.')
end

end


function [right_ascension,declination]=read_fallback_coordinates(filename)

header_file_name=resolve_header_file_name(filename,pwd);
header=get_header(header_file_name);
required_headers=get_required_headers_from_header(header);

if isfield(required_headers,'RA_2000') && ...
        isfield(required_headers,'DEC_2000')
    right_ascension=required_headers.RA_2000;
    declination=required_headers.DEC_2000;
elseif isfield(required_headers,'RA') && isfield(required_headers,'DEC')
    right_ascension=required_headers.RA;
    declination=required_headers.DEC;
else
    error(['Could not find RA/DEC or RA_2000/DEC_2000 in the FITS ' ...
        'header for barycentric correction.'])
end

end


function header_file_name=resolve_header_file_name(filename,starting_directory)

if isstring(filename), filename=char(filename); end
if isfile(filename)
    header_file_name=filename;
elseif isfile(fullfile(starting_directory,filename))
    header_file_name=fullfile(starting_directory,filename);
else
    header_file_name=filename;
end

end
