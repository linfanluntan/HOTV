% INPUTS
%------------------------------------
% file_list = {'20151005_0807431005A.nii';'20151005_0807431005AA.nii';'20151005_0807431005AB.nii';'20151005_0807431005AC.nii';'20151005_0807431005AD.nii';'20151005_0807431005AE.nii'};
% file_list = {'echo1.nii';'echo2.nii';'echo3.nii';'echo4.nii';'echo5.nii';'echo6.nii'};
file_list = {'10-2-e1-s.nii';'10-2-e2-s.nii';'10-2-e3-s.nii';'10-2-e4-s.nii';'10-2-e5-s.nii';'10-2-e6-s.nii'};
					% must point to valid nifti files
parameter_list = [1.05 2.2 3.35 4.5 5.65 6.8]';
					% units of ms or degrees
fit_type = 't2_exponential_plus_c'; 
					% options{'none','t2_linear_simple','t2_linear_weighted','t2_exponential','t2_linear_fast'
					%			't1_tr_fit','t1_fa_fit','t1_fa_linear_fit','t1_ti_exponential_fit'}
odd_echoes = 0;		% boolean, if selected only odd parameters will be
					% used for fit
rsquared_threshold = 0.2;
					% all fits with R^2 less than this set to -1
number_cpus = 4;	% not used if running on neuroecon
neuroecon = 0;		% boolean
output_basename = 'iron_nested';
					% base of output filename
data_order = 'xyzfile';% in what order is the data organized
					% options{'xynz','xyzn','xyzfile'}
tr = 20;			% units ms, only used for T1 FA fitting
email = '';
					% Email will be sent to this address on job completion
save_log       = 1;
email_log      = 0;
batch_log      = 0;
current_dir    = pwd;
log_name       = fullfile(current_dir,'log_t2.txt');
submit         = 1;
save_txt       = 1;
roi_list       = '';
fit_voxels     = 1;
xy_smooth_size = 0;


for m=size(file_list,1):-1:1
    testfile=cell2mat(file_list(m));
    file_list(m) = {fullfile(current_dir,testfile)};
end


cur_dataset.file_list = file_list;
cur_dataset.parameters = parameter_list;
parameter_list(:) = parameter_list;
cur_dataset.fit_type = fit_type;
cur_dataset.odd_echoes = odd_echoes;
cur_dataset.rsquared = rsquared_threshold;
cur_dataset.tr = tr;
cur_dataset.data_order = data_order;
cur_dataset.output_basename = output_basename;
cur_dataset.roi_list = roi_list;
cur_dataset.fit_voxels = fit_voxels;
cur_dataset.xy_smooth_size = xy_smooth_size;

JOB_struct(1).number_cpus = number_cpus ;
JOB_struct(1).neuroecon = neuroecon;
JOB_struct(1).email = email;
JOB_struct(1).batch_data = cur_dataset;
JOB_struct(1).save_log = save_log;
JOB_struct(1).email_log = email_log;
JOB_struct(1).batch_log = batch_log;
JOB_struct(1).current_dir = current_dir;
JOB_struct(1).log_name = log_name;
JOB_struct(1).submit = submit;
JOB_struct(1).save_txt = save_txt;
%------------------------------------

% Call Mapping Function
calculateMap(JOB_struct);
