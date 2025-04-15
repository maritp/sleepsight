function out = sleeplay_paths(local)

switch local
    case 1
        basedir = fullfile(filesep, 'Users', 'petzka', 'Documents', ...
            'projects', 'MPI', 'sleepsight_upload');
    case 0
        basedir = fullfile(filesep, 'Volumes', 'MPRG-Neurocode', ...
            'Data', 'sleepsight_2022_anikamarit');
end


%out.eegrawdir = fullfile(basedir, 'eeg', 'raw', 'all');
out.hypnodir = fullfile(basedir, 'eeg', 'sleepscoring', 'scored_rater');
out.eegdir = fullfile(basedir, 'eeg');
out.basedir = basedir;

% directory to fieldtrip 
ftv = 'fieldtrip-20220617';
ftdir = '/Users/petzka/Documents/projects/MPI';
out.fieldtripdir = fullfile(ftdir, ftv);
if contains(path,'fieldtrip')
    ft_defaults
else
    cd(out.fieldtripdir)
    ft_defaults
end


end
