

clc; clear all; close all;
%% Load phase data and mask:
% load('wrapphase_gre_3D.mat');
% load('mask_gre_3D.mat');
% phase_wrap = double(phase_wrap)*(2*pi)/256;
% % (for file size limitations on file exchange, the phase was stored as uint8)
% voxelSize = [.6 .6 .6];
% B0 = 3;
% gm = 2*pi*42.58e6;
% TE = 8.1e-3;
%
% % Unwrap the phase data:
% % phaseUnwrapped = laplacianUnwrap(phase_wrap, mask);
% % phaseUnwrapped = laplacianUnwrap(phase_wrap, mask);
%
%
% param.y=double(phase_wrap);
% recon_cs=param.y;


[namein, pathin, filterindex] = uigetfile({  '*.nii','NIFTI image (*.nii)'}, 'Select one or several files (using +CTRL or +SHIFT)','MultiSelect', 'on');

        if (iscell(namein))
            nbfile = size(namein,2);
        else
            nbfile = 1;
        end


img1=[];

for f = 1 : nbfile
    
    if(nbfile>1)
        filenamein = namein{f};
    else
        filenamein = namein;
    end
    
   disp(['Input file : ', fullfile(pathin, filenamein)])
   [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
    %    nout=[name_s 'updated' ext];
    %    pathout = pathin;
    %    disp(['Output file: ', fullfile(pathout, nout)])
    cname = fullfile(pathin, filenamein);
    VI=spm_vol(cname);
    
    
    if nbfile==1 && size(VI,1)>1
%         s=size(ima);
%         for vsn=1:size(VI,1)
            img1=spm_read_vols(VI);
            s=size(img1);
            sizedim=VI(1).dim;
        slicenum=sizedim(3);
            continue;
%         end
        
    else
        
        ima=spm_read_vols(VI);
        %    junk01=size(VI);
        %    gradients=junk01(1);
        sizedim=VI.dim;
        slicenum=sizedim(3);
        
        s=size(ima);
        ima=reshape(ima,s(1),s(2),s(3));
        %    imb =zeros(s(2),s(1),s(3));
        img1=cat(4,img1,ima);
    end
    
end



% for f = 1 : nbfile
%    
%    if(nbfile>1)
%        filenamein = namein{f};
%    else
%        filenamein = namein;
%    end
%    
%    disp(['Input file : ', fullfile(pathin, filenamein)])
%    [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
% %    nout=[name_s 'updated' ext];
% %    pathout = pathin;
% %    disp(['Output file: ', fullfile(pathout, nout)])  
%    cname = fullfile(pathin, filenamein);
%    VI=spm_vol(cname);
%    ima=spm_read_vols(VI);
% %    junk01=size(VI);
% %    gradients=junk01(1);
%    sizedim=VI.dim;
%    slicenum=sizedim(3);
%     
%    s=size(ima);
%    ima=reshape(ima,s(1),s(2),s(3));
% %    imb =zeros(s(2),s(1),s(3));
% img1=cat(4,img1,ima);
% 
% end




% nii=load_untouch_nii('C:\TACC\rock\30x560_1680dce_t1w.nii');


param=repmat({struct('y',[],'TV_dim1',[],'TV2_dim1',[],'TVWeight_dim1',[],'TVOP3D2D',[],'TVOPWeight',[],'TV2OP3D2D',[],'TV2OPWeight',[],'nite',[])}, 1,s(3));
parfor slicen=1:s(3) %1 %30
    param{slicen}.y=squeeze(double(img1(:,:,slicen,:)));
    recon_cs=param{slicen}.y;
    
    
    
    param{slicen}.TV_dim1=TV_Temp;  % temporal TV   2017-4-11
    param{slicen}.TVWeight_dim1=max(abs(recon_cs(:)))*0.03/2;%0.025;%.1;%0.25/2; % 0.03; temporal TV weight   2017-4-11
    
    param{slicen}.TV2_dim1=TV2_Temp;  % temporal TV   2017-4-11
    param{slicen}.TV2Weight_dim1=max(abs(recon_cs(:)))*0.03/2;%0.025;%.1;%0.25/2; % 0.03; temporal TV weight   2017-4-11
    
    param{slicen}.TVOP3D2D=TVOP3D2D; % spatial TV   2017-4-11
    param{slicen}.TVOPWeight=max(abs(recon_cs(:)))*0.7*.001/2;%param.TVWeight_dim1*0.09*0.09/0.025;%.1;%0.25/2; % 0.03  % spatial TV weight   2017-4-11
    
    param{slicen}.TV2OP3D2D=TV2OP3D2D; % spatial TV2   2017-4-11
    param{slicen}.TV2OPWeight=max(abs(recon_cs(:)))*0.3*.001/2;%param.TVOPWeight;%param.TVOPWeight*0.09;%.1;%0.25/2; % 0.03 % spatial TV2 weight   2017-4-11
    
    param{slicen}.nite = 50; %70; %100; %150; %300; %600;%40;%23;%100;% 35;%15;%15;%35;% 15, 25 iteration number   2017-4-11
    
    
    % param.y=double(phase_wrap);
    % recon_cs=param.y;
    
    clc
    tic
    recon_css = TPTVTV2TV2(recon_cs,param{slicen});
    
    % recon_csss = TPTVTV2SVT(recon_cs,param);
    
    [sx,sy,sz] = size(recon_css);
    img1(:,:,slicen,:)=reshape(int16(recon_css),sx,sy,1,sz);
end
clear param;



for f = 1 : nbfile
   
   if(nbfile>1)
       filenamein = namein{f};
   else
       filenamein = namein;
   end
   
   disp(['Input file : ', fullfile(pathin, filenamein)])
   [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
   
   
   nii=load_untouch_nii(fullfile(pathin, filenamein));
   nii.img=squeeze(img1(:,:,:,f));
   save_untouch_nii(nii,fullfile(pathin, ['TV2flt' filenamein]));
   
% %    nout=[name_s 'updated' ext];
% %    pathout = pathin;
% %    disp(['Output file: ', fullfile(pathout, nout)])  
%    cname = fullfile(pathin, filenamein);
%    VI=spm_vol(cname);
%    ima=spm_read_vols(VI);
% %    junk01=size(VI);
% %    gradients=junk01(1);
%    sizedim=VI.dim;
%    slicenum=sizedim(3);
%     
%    s=size(ima);
%    ima=reshape(ima,s(1),s(2),s(3));
% %    imb =zeros(s(2),s(1),s(3));
% img1=cat(4,img1,ima);

end

% save_untouch_nii(nii,'C:\TACC\rock\30x560_1680dce_t1wf001.nii')

