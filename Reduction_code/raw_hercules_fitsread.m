function data=raw_hercules_fitsread(filename)

%Read FITS header
header=get_header(filename);
headers=get_required_headers_from_header(header);

%Validate header geometry
validate_raw_hercules_header(headers,filename);

%Open FITS file
file_identifier=fopen(filename,'r','ieee-be');
if file_identifier==-1, error(['Could not open FITS file: ' filename]); end
cleanup_file=onCleanup(@() fclose(file_identifier));

%Find start of FITS image data
header_byte_count=find_fits_data_start_byte(file_identifier,filename);

%Read image data
fseek_status=fseek(file_identifier,header_byte_count,'bof');
if fseek_status~=0, error(['Could not seek to image data in FITS file: ' filename]); end
number_image_pixels=headers.NAXIS1*headers.NAXIS2;
raw_data=fread(file_identifier,number_image_pixels,'int16=>double');
if numel(raw_data)~=number_image_pixels, error(['File could not be read completely: ' filename]); end

%make reading offset disappear
raw_data=raw_data-min(raw_data);

%Reshape image
try
    data=reshape(raw_data,[headers.NAXIS1 headers.NAXIS2 1]);
catch reshape_error
    error(['file could not be read: ' filename '. ' reshape_error.message])
end

%trimming poor areas of the image
if size(data,1)<35 || size(data,2)<60, error(['Image is too small to trim bad detector regions: ' filename]); end
data(1:34,:)=[];   %top and bottom bad parts (read in to be all at the start)
data(:,end-29:end)=[];  %edges
data(:,1:30)=[];

%rotate image so that it is blue to red bottom to top and left to right
data=fliplr(data);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_raw_hercules_header(headers,filename)

%Check image dimensions
if ~isfield(headers,'NAXIS1') || isempty(headers.NAXIS1) || ~isfinite(headers.NAXIS1) || headers.NAXIS1<=0, error(['Missing or invalid NAXIS1 in FITS file: ' filename]); end
if ~isfield(headers,'NAXIS2') || isempty(headers.NAXIS2) || ~isfinite(headers.NAXIS2) || headers.NAXIS2<=0, error(['Missing or invalid NAXIS2 in FITS file: ' filename]); end
if headers.NAXIS1~=round(headers.NAXIS1) || headers.NAXIS2~=round(headers.NAXIS2), error(['NAXIS1 and NAXIS2 must be integer values in FITS file: ' filename]); end

end

function header_byte_count=find_fits_data_start_byte(file_identifier,filename)

%Find FITS END card
fseek(file_identifier,0,'bof');
header_byte_count=0;
end_found=false;
while ~end_found
    header_block=fread(file_identifier,2880,'*char')';
    if isempty(header_block), error(['Reached end of file before FITS END card was found: ' filename]); end
    if numel(header_block)~=2880, error(['Incomplete FITS header block in file: ' filename]); end
    header_byte_count=header_byte_count+2880;
    for card_index=1:36
        card_start=(card_index-1)*80+1;
        card_stop=card_index*80;
        header_card=char(header_block(card_start:card_stop));
        if startsWith(header_card,'END')
            end_found=true;
            break
        end
    end
end

end