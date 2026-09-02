function [star,identified]=find_barycentric_catalogue_target(target_name,catalogue_path)
%FIND_BARYCENTRIC_CATALOGUE_TARGET Resolve supported stellar names.
%   HD, HIP, HR, Bayer, and Flamsteed designations are supported, including
%   HRSP-style underscores and optional single-letter components. Examples:
%   HD_209295, HIP 108976, HR 5459, Alpha Centauri A, Alp_Cen_A, 9 CMa.
%   An unidentified or ambiguous target issues a warning and returns false.

if nargin<2 || isempty(catalogue_path)
    catalogue_path=fullfile(fileparts(mfilename('fullpath')), ...
        'barycentric_stellar_catalogue.mat');
end
if ~isfile(catalogue_path)
    error(['Converted HRSP catalogue not found: %s\nRun ' ...
        'convert_hrsp_catalogues_to_mat first.'],catalogue_path)
end

persistent cached_catalogue cached_path
if isempty(cached_catalogue) || isempty(cached_path) || cached_path~=string(catalogue_path)
    loaded=load(catalogue_path,'catalogue');
    cached_catalogue=loaded.catalogue;
    cached_path=string(catalogue_path);
end
catalogue=cached_catalogue;
if ~isfield(catalogue,'bright') || ~isfield(catalogue,'ccdm')
    error(['The converted catalogue predates expanded name lookup. Run ' ...
        'convert_hrsp_catalogues_to_mat again.'])
end

target_text=normalise_name(target_name);
[designation,parsed]=parse_designation(target_text,catalogue);
if ~parsed
    [star,identified]=not_identified(target_text, ...
        'name is not a recognised HD, HIP, HR, Bayer, or Flamsteed designation');
    return
end

if designation.kind=="USR"
    [star,identified]=not_identified(target_text, ...
        'USR entries are observation-specific and are not stored in the converted catalogue');
    return
end

if any(designation.kind==["HR","BAYER","FLAMSTEED"])
    [hd,component,reason]=bright_to_hd(designation,catalogue.bright);
    if ~isempty(reason)
        [star,identified]=not_identified(target_text,reason);
        return
    end
    designation.kind="HD";
    designation.number=hd;
    designation.component=component;
end

[record,row,reason]=resolve_hip_record(designation,catalogue);
if ~isempty(reason)
    [star,identified]=not_identified(target_text,reason);
    return
end

star=struct();
fields=fieldnames(record);
for field_index=1:numel(fields)
    value=record.(fields{field_index});
    star.(fields{field_index})=value(row,:);
end
star.reference_frame=catalogue.reference_frame;
star.equinox=catalogue.equinox;
star.epoch=catalogue.epoch;
star.input_designation=string(target_name);
identified=true;
end


function [designation,ok]=parse_designation(text,catalogue)

designation=struct('kind',"",'number',NaN,'component',"", ...
    'bayer',"",'sequence',0,'constellation',"");
ok=false;

tokens=regexp(text,'^(HD|HIP|HR|USR)\s*(\d+)(?:\s+([A-Z]))?$','tokens','once');
if ~isempty(tokens)
    designation.kind=string(tokens{1});
    designation.number=str2double(tokens{2});
    if numel(tokens)>=3, designation.component=string(tokens{3}); end
    ok=true;
    return
end

tokens=regexp(text,'^(\d+)\s+(.+)$','tokens','once');
if ~isempty(tokens)
    [constellation,component,matched]=parse_constellation(tokens{2},catalogue);
    if matched
        designation.kind="FLAMSTEED";
        designation.number=str2double(tokens{1});
        designation.constellation=constellation;
        designation.component=component;
        ok=true;
    end
    return
end

greek_full=upper(catalogue.greek.full);
greek_abbreviation=upper(catalogue.greek.abbreviation);
best_length=0;
greek_row=0;
for row=1:numel(greek_full)
    aliases=[greek_full(row),greek_abbreviation(row)];
    for alias=aliases
        length_alias=strlength(alias);
        if length_alias<=best_length || ~startsWith(text,alias), continue; end
        if strlength(text)>length_alias
            following=extractBetween(text,length_alias+1,length_alias+1);
            if ~isempty(regexp(char(following),'[A-Z]','once')), continue; end
        end
        best_length=length_alias;
        greek_row=row;
    end
end
if greek_row==0, return; end

remainder=strtrim(extractAfter(text,best_length));
sequence_token=regexp(remainder,'^(\d)(.*)$','tokens','once');
sequence=0;
if ~isempty(sequence_token)
    sequence=str2double(sequence_token{1});
    remainder=strtrim(string(sequence_token{2}));
end
[constellation,component,matched]=parse_constellation(remainder,catalogue);
if ~matched, return; end
designation.kind="BAYER";
designation.bayer=catalogue.greek.abbreviation(greek_row);
designation.sequence=sequence;
designation.constellation=constellation;
designation.component=component;
ok=true;
end


function [abbreviation,component,matched]=parse_constellation(text,catalogue)

text=string(text);
genitive=upper(catalogue.constellations.genitive);
short=upper(catalogue.constellations.abbreviation);
best_length=0;
best_row=0;
for row=1:numel(genitive)
    aliases=[genitive(row),short(row)];
    for alias=aliases
        length_alias=strlength(alias);
        if length_alias<=best_length || ~startsWith(text,alias), continue; end
        if strlength(text)>length_alias && extractBetween(text,length_alias+1,length_alias+1)~=" "
            continue
        end
        remainder=strtrim(extractAfter(text,length_alias));
        if strlength(remainder)>1, continue; end
        if strlength(remainder)==1 && isempty(regexp(char(remainder),'^[A-Z]$','once'))
            continue
        end
        best_length=length_alias;
        best_row=row;
    end
end
matched=best_row>0;
if ~matched
    abbreviation="";
    component="";
    return
end
abbreviation=catalogue.constellations.abbreviation(best_row);
component=strtrim(extractAfter(text,best_length));
end


function [hd,component,reason]=bright_to_hd(designation,bright)

component=designation.component;
reason='';
switch designation.kind
    case "HR"
        if designation.number<1 || designation.number>numel(bright.hr)
            hd=NaN; reason='HR number is outside the Bright Star Catalogue'; return
        end
        row=designation.number;
    case "BAYER"
        matches=bright.valid & strcmpi(bright.bayer,designation.bayer) & ...
            strcmpi(bright.constellation,designation.constellation);
        if designation.sequence~=0, matches=matches & bright.sequence==designation.sequence; end
        if strlength(component)>0, matches=matches & contains(bright.component,component); end
        rows=find(matches);
        if numel(rows)~=1
            hd=NaN; reason=match_reason(rows,'Bright Star Catalogue'); return
        end
        row=rows(1);
    case "FLAMSTEED"
        matches=bright.valid & bright.flamsteed==designation.number & ...
            strcmpi(bright.constellation,designation.constellation);
        if strlength(component)>0, matches=matches & contains(bright.component,component); end
        rows=find(matches);
        if numel(rows)~=1
            hd=NaN; reason=match_reason(rows,'Bright Star Catalogue'); return
        end
        row=rows(1);
    otherwise
        error('Unexpected Bright Star designation kind.')
end
hd=double(bright.hd(row));
if hd<1, reason='Bright Star Catalogue entry has no HD mapping'; end
end


function [record,row,reason]=resolve_hip_record(designation,catalogue)

record=struct(); row=[]; reason='';
component=designation.component;

if designation.kind=="HD"
    ccdm_matches=catalogue.ccdm.hd>0 & ...
        (catalogue.ccdm.hd==designation.number | ...
        (catalogue.ccdm.double_hd & catalogue.ccdm.hd==designation.number+1));
    ccdm_rows=find(ccdm_matches);
    if numel(ccdm_rows)>1
        reason='multiple matching CCDM records'; return
    elseif numel(ccdm_rows)==1
        ccdm_row=ccdm_rows(1);
        ccdm_component=char(catalogue.ccdm.component(ccdm_row));
        if numel(ccdm_component)>=2
            expected_component=string(ccdm_component(2));
            if strlength(component)>0 && component~=expected_component
                reason='component is inconsistent with the CCDM mapping'; return
            end
            matches=catalogue.components.ccdm==catalogue.ccdm.ccdm(ccdm_row) & ...
                startsWith(catalogue.components.component,expected_component);
            rows=find(matches);
            if numel(rows)~=1
                reason=match_reason(rows,'Hipparcos component catalogue'); return
            end
            record=catalogue.components; row=rows(1); return
        end
    end
    matches=catalogue.main.valid & catalogue.main.hd==designation.number;
else
    matches=catalogue.main.valid & catalogue.main.hip==designation.number;
end

rows=find(matches);
if numel(rows)~=1
    reason=match_reason(rows,'Hipparcos main catalogue'); return
end
main_row=rows(1);
hip=catalogue.main.hip(main_row);

available_components=unique([catalogue.main.component(catalogue.main.hip==hip); ...
    catalogue.components.component(catalogue.components.hip==hip)]);
available_components=join(available_components,'');
if strlength(component)==0
    if catalogue.main.annex(main_row)=="C" && catalogue.main.number_components(main_row)>1
        reason='multiple Hipparcos components; a component letter is required'; return
    end
else
    if ~contains(available_components,component)
        reason='component is not present in the Hipparcos catalogue'; return
    end
    if catalogue.main.annex(main_row)=="C"
        component_rows=find(catalogue.components.hip==hip & ...
            catalogue.components.component==component);
        if numel(component_rows)~=1
            reason=match_reason(component_rows,'Hipparcos component catalogue'); return
        end
        record=catalogue.components; row=component_rows(1); return
    end
end
record=catalogue.main;
row=main_row;
end


function reason=match_reason(rows,catalogue_name)

if isempty(rows)
    reason=['no matching valid record in the ' catalogue_name];
else
    reason=['multiple matching records in the ' catalogue_name];
end
end


function text=normalise_name(value)

text=upper(strtrim(string(value)));
text=regexprep(text,'[^A-Z0-9]+',' ');
text=strtrim(regexprep(text,'\s+',' '));
end


function [star,identified]=not_identified(target_text,reason)

star=[];
identified=false;
warning('MEGARA:BarycentricTargetNotIdentified', ...
    'Target "%s" was not identified: %s.',target_text,reason);
end
