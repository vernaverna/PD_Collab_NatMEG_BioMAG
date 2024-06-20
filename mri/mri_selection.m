clear all; clc

fold = ''  %provide data folder path
protos = dir(fold);
protos = protos([protos.isdir]);

proto_names = {};
proto_folders = {};
noDCMfiles = [];
for f = 3:length(protos)
    files = dir([fold filesep protos(f).name filesep '*.dcm']);

    if ~isempty(files)
        dcminfo = dicominfo([fold filesep protos(f).name filesep files(1).name]);
        if isfield(dcminfo, 'ProtocolName')
            proto_n = dcminfo.ProtocolName;
            proto_n = proto_n(find(~isspace(proto_n))); %remove spaces
            proto_names{end+1} = proto_n;
            proto_folders{end+1} = protos(f).name;
            noDCMfiles = [noDCMfiles length(files)];        

        end
    end
end

proto_names'
noDCMfiles'
proto_folders'

%%

in_file = [fold filesep proto_folders{9} filesep '1.1.1.1.11111111.111111111111111.dcm']; %provide one dicome file from the selected sequence folder
out_file = ['/projects/sub-111/sub-111_T1w.nii.gz']; %provide output folder and file name
mri_convert_command = ['mri_convert ' in_file ' ' out_file];
[status,cmdout] = system(mri_convert_command);
disp(cmdout)
