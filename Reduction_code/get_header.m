function header=get_header(filename)

%Validate input
if isstring(filename), filename=char(filename); end
if ~ischar(filename) || isempty(filename), error('filename must be a non-empty character vector or string scalar.'); end
if ~isfile(filename), error(['Header file not found: ' filename]); end

%Open FITS file
file_identifier=fopen(filename,'r','ieee-be');
if file_identifier==-1, error(['Could not open file: ' filename]); end
cleanup_file=onCleanup(@() fclose(file_identifier));

%Read FITS header cards
header={};
card_index=0;
end_found=false;
while ~end_found
    header_block=fread(file_identifier,2880,'*char')';
    if isempty(header_block), error(['Reached end of file before FITS END card was found: ' filename]); end
    if numel(header_block)~=2880, error(['Incomplete FITS header block in file: ' filename]); end
    for block_card_index=1:36
        card_index=card_index+1;
        card_start=(block_card_index-1)*80+1;
        card_stop=block_card_index*80;
        header{card_index}=char(header_block(card_start:card_stop));
        if startsWith(header{card_index},'END')
            end_found=true;
            break
        end
    end
end

end

