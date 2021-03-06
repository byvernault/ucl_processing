%% Matlab Code for AutoSpider for ADC_MAP
% add path where matlab code is from inputs:
addpath(genpath('${matlab_code}'));

% Paths:
mkdir('${temp_dir}/outputs');
mkdir('${temp_dir}/OsiriX');

% Get dicoms:
dicoms = '${dicoms}';
[dcmpath, dcmname, dcmext] = fileparts(dicoms);
if strcmp(dcmext, '.zip')
	cmd = sprintf('cd %s && unzip ${dicoms} > /dev/null', dcmpath)
	system(cmd, '-echo');
end

% Call function:
XNAT_ADCsMaps(dcmpath, '${temp_dir}/outputs', '${temp_dir}/adc_map_report.pdf');

% copy the outputs files:
files = dir('${temp_dir}/outputs/ADC/*.dcm');
for file = files'
    copyfile(['${temp_dir}/outputs/ADC/' file.name], '${temp_dir}/OsiriX')
end

files = dir('${temp_dir}/outputs/ADCNoZero/*.dcm');
for file = files'
    copyfile(['${temp_dir}/outputs/ADCNoZero/' file.name], '${temp_dir}/OsiriX')
end
files = dir('${temp_dir}/outputs/ADCfast/*.dcm');
for file = files'
    copyfile(['${temp_dir}/outputs/ADCfast/' file.name], '${temp_dir}/OsiriX')
end
files = dir('${temp_dir}/outputs/ADCslow/*.dcm');
for file = files'
    copyfile(['${temp_dir}/outputs/ADCslow/' file.name], '${temp_dir}/OsiriX')
end

% zip the resource OsiriX folder:
cmd = 'cd ${temp_dir}/OsiriX && zip -r osirix.zip * > /dev/null && mv osirix.zip .. && cd ..'
system(cmd, '-echo');