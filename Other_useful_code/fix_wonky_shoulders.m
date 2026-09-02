%FIX_WONKY_SHOULDERS Manually renormalise a legacy coarse CCF baseline.
% This interactive maintenance utility is retained for old CCF products;
% current regularised-LSD reductions do not normally need it. Set the target
% below, run beside its final_data and ccf files, and click baseline points
% when prompted. Corrected CCF and FAMIAS files are written to a new folder.

clearvars
objectname='hd_175362';


load ("final_data_" + objectname + ".mat")
load ("ccf_" + objectname + ".mat")
fullwave=finaldata.fullwave;
fullint=finaldata.fullint;
%jd=finaldata.jd;


figure
plot(coarse_ccf.velocity,coarse_ccf.intensity)
[x,y] = getpts;
p = polyfit(x,y,1);
f1 = polyval(p,coarse_ccf.velocity);
norm_int=coarse_ccf.intensity./f1;
hold on
plot(coarse_ccf.velocity,norm_int)
coarse_ccf.intensity=norm_int;


%getting weights from a noise estimate and producing a summed spectrum

%Find Nans and set to 1
G=isnan(fullint);
fullint(G)=1;

intfilt=medfilt1(fullint(:,1:10000)',11)';
resid=fullint(:,1:10000)-intfilt;
weight=1./std(resid,0,2);

U=coarse_ccf.velocity;
LSDy=coarse_ccf.intensity;
jd=coarse_ccf.jd;


output_directory=fullfile('Reduced_Data',objectname, ...
    ['Wonky_fix_CCF_FAMIAS_files_' objectname(4:end)]);
if ~isfolder(output_directory)
    mkdir(output_directory)
end
save(fullfile(output_directory,['ccf_' objectname '_wonky_fix.mat']), ...
    'coarse_ccf')
write_files_for_FAMIAS(jd,U,LSDy,weight, ...
    fullfile(output_directory,['times_' objectname(4:end) '.txt']), ...
    fullfile(output_directory,[objectname(4:end) '_']))
