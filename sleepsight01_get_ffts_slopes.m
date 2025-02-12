% SLEEPSIGHT. 
%
%
% obtain ffts & slopes
%
%
% c Marit Petzka, University of Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%%
local_ = 1; % is data local or retrieved from server 

%% getting directories
gitdir = fullfile(filesep, 'Users', 'petzka', ...
    'Documents', 'GitHub', 'sleepsight');
addpath(genpath(gitdir))

dirs = sleepsight_paths(local_);

%% where relevant data is
datfolder = 'data_ready4analyses';
datdir = fullfile(dirs.eegdir, datfolder);
files_ = dir(fullfile(datdir, '*.mat'));

%% where output should be saved
savename = 'power_spectra_1_45hz';
savefolder = 'fft'; % data
savedir = dirs.eegdir;

%% settings
a = {files_.name};
settings_.nsub = cellfun(@(x) x(1:4), a, 'UniformOutput',false);
settings_.length = 6; % length in seconds of 'trials' for FFT
settings_.segoverlap = .5; % default = 0. number between 0 and 1
settings_.FOI = [1 45];
settings_.freqres = 0.2; %frequency resolution
settings_.de = 1; % doing detrend and demeaning on whole data
% settings_.fft_slope = 1; %0 = no. 1 = fooof. 
settings_.rej_muc = 1; % reject artifacts based on visuals
settings_.scoring = cell(1,numel(settings_.nsub));

%% predefine
fft_ = cell(1,numel(settings_.nsub));

fft_slope.osci = cell(1,numel(settings_.nsub));
fft_slope.mixd = cell(1,numel(settings_.nsub));
fft_slope.frac = cell(1,numel(settings_.nsub));

%% lets gooooo
tic

for isub = 1:numel(settings_.nsub)

    fileOI = files_(isub).name;

    fprintf('----------subject %s. %d/%d\n', fileOI(1:4), isub, numel(settings_.nsub))

    eeg = load(fullfile(datdir, fileOI));
    dat = eeg.dat;
    artifacts = dat.artifacts;

    %% check scoring.
    settings_.scoring{isub} = eeg.dat.scoring{1};

    %% detrend & demean data 
    if settings_.de == 1
        cfg = [];
        cfg.detrend = 'yes';

        dat = ft_preprocessing(cfg, dat);
        dat.artifacts = artifacts;

    end
    %% --- cut into x sec trials
    cfg = [];
    cfg.length = settings_.length;
    cfg.overlap = settings_.segoverlap;
    dat_trl = ft_redefinetrial(cfg, dat);

    %--- give all trials the same time axis
    for itime = 1:numel(dat_trl.time)
        dat_trl.time{itime} = 1/dat_trl.fsample:1/dat_trl.fsample:settings_.length;
    end

    %% exclude trials containing artifacts

    if settings_.rej_muc == 1

        if ~isempty(dat.artifacts)
            art_def = hvn_createBnrySignal(dat.artifacts,...
                dat_trl.sampleinfo(end,2)); % get logical vector

            idx_trl_excl = zeros(1,size(dat_trl.sampleinfo,1));

            for itrl = 1:size(dat_trl.sampleinfo,1)
                a = art_def(dat_trl.sampleinfo(itrl,1):dat_trl.sampleinfo(itrl,2));
                if sum(a) > 0
                    idx_trl_excl(itrl) = 1;
                end
            end

            dat_trl_tmp = dat_trl;
            dat_trl_tmp.trial = dat_trl_tmp.trial(~idx_trl_excl);
            dat_trl_tmp.time = dat_trl_tmp.time(~idx_trl_excl);
            dat_trl_tmp.sampleinfo = dat_trl_tmp.sampleinfo(~idx_trl_excl',:);
            dat_trl = dat_trl_tmp;
        end
    end

    %% fft
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.foi = settings_.FOI(1):settings_.freqres:settings_.FOI(2);
    %cfg.foilim = settings_.FOI;
    %cfg.foilim = settings_.FOI;
    cfg.taper = 'hanning';
    cfg.pad = 'nextpow2';
    cfg.keeptrials = 'yes';

    fft_tmp = ft_freqanalysis(cfg, dat_trl);

    %% save
    %fft_tmp = rmfield(fft_tmp, {'powspctrm' 'cumsumcnt' 'cumtapcnt'});
    fft_{isub} = fft_tmp;

    %% FOOOOOOOFFFF

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fooof';
    cfg.foi = settings_.FOI(1):settings_.freqres:settings_.FOI(2);
    cfg.taper = 'hanning';
    cfg.pad = 'nextpow2';

    fft_slope.mixd{isub} = ft_freqanalysis(cfg, dat_trl);

    cfg.output = 'fooof_peaks';
    fft_slope.osci{isub} = ft_freqanalysis(cfg, dat_trl);

    cfg.output = 'fooof_aperiodic';
    fft_slope.frac{isub} = ft_freqanalysis(cfg, dat_trl);


end

settings_.fsample = dat.fsample;

%% tic toc
fprintf('---------- overall time needed: %d\n', toc)

%% save
if ~exist(fullfile(savedir, savefolder))
    mkdir(fullfile(savedir, savefolder))
end

save(fullfile(savedir, savefolder, savename),...
    'fft_', 'fft_slope', 'settings_', '-v7.3')
