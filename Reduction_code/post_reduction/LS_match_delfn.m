function [newdelfn,thefits,waves]=LS_match_delfn(wavelength,spec,delfn,vsini,intrinFWHM,epsx)

%Validate inputs
validate_LS_match_delfn_inputs(wavelength,spec,delfn,vsini,intrinFWHM,epsx);

%Prepare arrays
wavelength=wavelength(:);
spec=spec(:);
delfn=sortrows(delfn,1);
newdelfn=[];
thefits=struct();
waves=struct();

%Set fitting regions
start=delfn(1,1);
cnt=0;
crossover=7; %size of crossover wavelength in each piece fitted (7 is good enough)
stepsize=50; %size of wavelength region to be fitted at one time (50 is pretty good with reasonable fitting times)

%Fit line strengths region by region
while start<delfn(end,1)
    finish=start+stepsize;
    delind=find(delfn(:,1)>=start & delfn(:,1)<finish);
    if isempty(delind)
        %no lines in this region
        start=finish;
        continue
    end
    waveind=find(wavelength>start-1 & wavelength<finish+1);
    if isempty(waveind)
        start=finish-crossover;
        continue
    end
    cnt=cnt+1;
    tdelfn=delfn(delind,:);
    wave=wavelength(waveind);
    sp=spec(waveind);
    [final_delfns,residual,thefit]=leastsquaresfit_spec_w_broad_delfns(wave,sp,tdelfn,vsini,intrinFWHM,epsx);
    thefits.(['f' num2str(cnt)])=thefit;
    waves.(['w' num2str(cnt)])=wave;
    if cnt==1
        %chop end (not start) off since there is overlap
        wantedinds=find(final_delfns(:,1)<finish-(crossover/2));
        final_delfns=final_delfns(wantedinds,:);
        newdelfn=final_delfns;
    else
        %chop both ends off to get rid of poor end fitting (this is why there is
        %overlap between fits)
        wantedinds=find(final_delfns(:,1)>=start+(crossover/2) & final_delfns(:,1)<finish-(crossover/2));
        final_delfns=final_delfns(wantedinds,:);
        newdelfn=cat(1,newdelfn,final_delfns);
    end
    start=finish-crossover;
end

%Check output
if isempty(newdelfn), warning('LS_match_delfn did not produce any matched delta-function lines. Returning input delfn.'); newdelfn=delfn; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_LS_match_delfn_inputs(wavelength,spec,delfn,vsini,intrinFWHM,epsx)

%Check spectrum inputs
if ~isnumeric(wavelength) || ~isvector(wavelength), error('wavelength must be a numeric vector.'); end
if ~isnumeric(spec) || ~isvector(spec), error('spec must be a numeric vector.'); end
if numel(wavelength)~=numel(spec), error('wavelength and spec must have the same number of elements.'); end
if any(~isfinite(wavelength(:))) || any(~isfinite(spec(:))), error('wavelength and spec must contain only finite values.'); end
if any(diff(wavelength(:))<=0), error('wavelength must be strictly increasing.'); end

%Check delta-function input
if ~isnumeric(delfn) || size(delfn,2)~=2, error('delfn must be a numeric two-column matrix: wavelength and weight.'); end
if isempty(delfn), error('delfn is empty.'); end
if any(~isfinite(delfn(:))), error('delfn must contain only finite values.'); end

%Check broadening inputs
if ~isnumeric(vsini) || ~isscalar(vsini) || ~isfinite(vsini) || vsini<0, error('vsini must be a finite scalar value >= 0.'); end
if ~isnumeric(intrinFWHM) || ~isscalar(intrinFWHM) || ~isfinite(intrinFWHM) || intrinFWHM<0, error('intrinFWHM must be a finite scalar value >= 0.'); end
if ~isnumeric(epsx) || ~isscalar(epsx) || ~isfinite(epsx) || epsx<0 || epsx>1, error('epsx must be a finite scalar value between 0 and 1.'); end

end