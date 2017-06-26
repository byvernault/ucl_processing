%% Matlab Code for AutoSpider for Sample GM Segment using SPM
% add path where matlab code is from inputs:
addpath(genpath('${matlab_code}'));

% Get dicoms:
nii_image = '${t1}';
[niipath, niiname, niiext] = fileparts(nii_image);
if strcmp(niiext, '.gz')
	cmd = sprintf('gzip -d ${t1}');
	system(cmd, '-echo');
end

% Call function:
sample_GM_segment([niipath '/' niiname], '${spm}');

% Generate a report PDF for display on XNAT:
fl = dir([niipath '/c1*.nii']);
c1file = [niipath '/' fl(1).name];
fl = dir([niipath '/c2*.nii']);
c2file = [niipath '/' fl(1).name];
fl = dir([niipath '/c3*.nii']);
c3file = [niipath '/' fl(1).name];
fl = dir([niipath '/m*.nii']);
mfile = [niipath '/' fl(1).name];

% open figure:
fig = figure(1);

% load the images:
nii = load_nii(c1file);
img = nii.img();
ind = 0;

% Size of the data
img_size = size(img);
s1 = round(img_size(3)*2/6);
s2 = round(img_size(3)*3/8);
s3 = round(img_size(3)*5/8);
s4 = round(img_size(3)*4/6);

subplot(4,4,4*ind+1);
min_v = min(min(img(:,:,s1)));
max_v = max(max(img(:,:,s1)));
imshow(img(:,:,s1), [min_v, max_v]);
title('C1');
xlabel(['Slice: ' num2str(s1)]);

subplot(4,4,4*ind+2);
min_v = min(min(img(:,:,s2)));
max_v = max(max(img(:,:,s2)));
imshow(img(:,:,s2), [min_v, max_v]);
xlabel(['Slice: ' num2str(s2)]);

subplot(4,4,4*ind+3);
min_v = min(min(img(:,:,s3)));
max_v = max(max(img(:,:,s3)));
imshow(img(:,:,s3), [min_v, max_v]);
xlabel(['Slice: ' num2str(s3)]);

subplot(4,4,4*ind+4);
min_v = min(min(img(:,:,s4)));
max_v = max(max(img(:,:,s4)));
imshow(img(:,:,s4), [min_v, max_v]);
xlabel(['Slice: ' num2str(s4)]);

% load the images:
nii = load_nii(c2file);
img = nii.img();
ind = 1;

% Size of the data
img_size = size(img);
s1 = round(img_size(3)*2/6);
s2 = round(img_size(3)*3/8);
s3 = round(img_size(3)*5/8);
s4 = round(img_size(3)*4/6);

subplot(4,4,4*ind+1);
min_v = min(min(img(:,:,s1)));
max_v = max(max(img(:,:,s1)));
imshow(img(:,:,s1), [min_v, max_v]);
title('C2');
xlabel(['Slice: ' num2str(s1)]);

subplot(4,4,4*ind+2);
min_v = min(min(img(:,:,s2)));
max_v = max(max(img(:,:,s2)));
imshow(img(:,:,s2), [min_v, max_v]);
xlabel(['Slice: ' num2str(s2)]);

subplot(4,4,4*ind+3);
min_v = min(min(img(:,:,s3)));
max_v = max(max(img(:,:,s3)));
imshow(img(:,:,s3), [min_v, max_v]);
xlabel(['Slice: ' num2str(s3)]);

subplot(4,4,4*ind+4);
min_v = min(min(img(:,:,s4)));
max_v = max(max(img(:,:,s4)));
imshow(img(:,:,s4), [min_v, max_v]);
xlabel(['Slice: ' num2str(s4)]);

% load the images:
nii = load_nii(c3file);
img = nii.img();
ind = 2;

% Size of the data
img_size = size(img);
s1 = round(img_size(3)*2/6);
s2 = round(img_size(3)*3/8);
s3 = round(img_size(3)*5/8);
s4 = round(img_size(3)*4/6);

subplot(4,4,4*ind+1);
min_v = min(min(img(:,:,s1)));
max_v = max(max(img(:,:,s1)));
imshow(img(:,:,s1), [min_v, max_v]);
title('C3');
xlabel(['Slice: ' num2str(s1)]);

subplot(4,4,4*ind+2);
min_v = min(min(img(:,:,s2)));
max_v = max(max(img(:,:,s2)));
imshow(img(:,:,s2), [min_v, max_v]);
xlabel(['Slice: ' num2str(s2)]);

subplot(4,4,4*ind+3);
min_v = min(min(img(:,:,s3)));
max_v = max(max(img(:,:,s3)));
imshow(img(:,:,s3), [min_v, max_v]);
xlabel(['Slice: ' num2str(s3)]);

subplot(4,4,4*ind+4);
min_v = min(min(img(:,:,s4)));
max_v = max(max(img(:,:,s4)));
imshow(img(:,:,s4), [min_v, max_v]);
xlabel(['Slice: ' num2str(s4)]);

% load the images:
nii = load_nii(mfile);
img = nii.img();
ind = 3;

% Size of the data
img_size = size(img);
s1 = round(img_size(3)*2/6);
s2 = round(img_size(3)*3/8);
s3 = round(img_size(3)*5/8);
s4 = round(img_size(3)*4/6);

subplot(4,4,4*ind+1);
min_v = min(min(img(:,:,s1)));
max_v = max(max(img(:,:,s1)));
imshow(img(:,:,s1), [min_v, max_v]);
title('M');
xlabel(['Slice: ' num2str(s1)]);

subplot(4,4,4*ind+2);
min_v = min(min(img(:,:,s2)));
max_v = max(max(img(:,:,s2)));
imshow(img(:,:,s2), [min_v, max_v]);
xlabel(['Slice: ' num2str(s2)]);

subplot(4,4,4*ind+3);
min_v = min(min(img(:,:,s3)));
max_v = max(max(img(:,:,s3)));
imshow(img(:,:,s3), [min_v, max_v]);
xlabel(['Slice: ' num2str(s3)]);

subplot(4,4,4*ind+4);
min_v = min(min(img(:,:,s4)));
max_v = max(max(img(:,:,s4)));
imshow(img(:,:,s4), [min_v, max_v]);
xlabel(['Slice: ' num2str(s4)]);

saveas(fig, '${temp_dir}/sample_gm_report.pdf');

% copy the outputs files:
files = dir([niipath '/c*.nii']);
for file = files'
    copyfile([niipath '/' file.name], '${temp_dir}');
end
% m file:
copyfile(mfile, '${temp_dir}');
files = dir([niipath '/*.mat']);
for file = files'
    copyfile([niipath '/' file.name], '${temp_dir}');
end
files = dir([niipath '/rc*.nii']);
for file = files'
    copyfile([niipath '/' file.name], '${temp_dir}');
end

% Gzip nii images:
files = dir(['${temp_dir}' '/*.nii']);
for file = files'
	cmd = sprintf(['gzip ${temp_dir}/' file.name]);
	system(cmd, '-echo');
end
