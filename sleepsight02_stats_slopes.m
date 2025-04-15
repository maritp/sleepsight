% SLEEPSIGHT. 
%
%
% doing stats on slopes
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
filename = 'power_spectra_1_45hz';

datfolder = 'fft'; % data
datdir = fullfile(dirs.eegdir, datfolder);

filename_gr = 'sleepsight_behav.csv';
gr = readtable(fullfile(dirs.eegdir, filename_gr));

%% where output should be saved
savename = 'power_spectra_1_45hz';
savefolder = 'eeg4stats'; % data
savedir = dirs.eegdir;

%% getting layout
load(fullfile(gitdir, 'mpib-60Ch-layout.mat')) %layout
load(fullfile(gitdir, 'mpib-60Ch-neighbours.mat')) %neighbours
load(fullfile(gitdir, 'mpib-60Ch-elec.mat')) %electrode position

%% load
load(fullfile(datdir, filename))

%% settings
settings_.selFFTs = 'frac'; %analysis with: 'fft' = uncorr. fieldtrip.
                           %               'osci' = corr for 1/f 
                           %               'mixd' = uncorr. 
                           %               'frac' = 1/f fit. 

settings_.n2_thres = 2; % threshold for including n2 participants; in minutes
%% select ffts
if strcmpi(settings_.selFFTs, 'fft')
    dat = fft_;
else
    dat = fft_slope.(settings_.selFFTs);
end

%% goooooooo

for isub = 1:numel(dat)

    %% put the slope into the dat struct

    b = [dat{isub}.fooofparams.aperiodic_params];
    d = reshape(b, [numel(b)/numel(dat{isub}.label) ...
        numel(dat{isub}.label)]);
    d = d'; % ch x [offset, exponent]

    dat{isub}.freq = 1;
    dat{isub}.powspctrm = d(:,2);

    %% getting correspoding group assignments
    sub_no = str2num(settings_.nsub{isub}(2:end));

    sub_idx = find(gr.EEG_id == sub_no);

    gr_vars.insight(isub) = gr.insight(sub_idx); % insight
    gr_vars.sleep(isub) = gr.sleep_scoring(sub_idx); % EEG 
    gr_vars.sleep_report(isub) = gr.subjective_scoring(sub_idx); % subjective 

    gr_vars.behav_id(isub) = gr.behav_id(sub_idx); % behav id 
    gr_vars.eeg_id(isub) = gr.EEG_id(sub_idx); % eeg id 

    %% check for length of n2 
    n2_dur_sp = sum(settings_.scoring{isub} == 2);

    if n2_dur_sp >= settings_.n2_thres*60*settings_.fsample
        gr_vars.n2_long(isub) = 1;
    else
        gr_vars.n2_long(isub) = 0;
    end

    if n2_dur_sp > 0
        gr_vars.n2_dur(isub) = n2_dur_sp/settings_.fsample/60;
    end

end

%% insight likelihood
group_label = {'wake', 'N1', 'N2'};

for i = 1:numel(group_label)
    sprintf('group %s: insight.yes: %d, no: %d', ...
    group_label{i},...
    numel(gr_vars.behav_id((gr_vars.sleep == i-1 & gr_vars.insight == 1))),...
    numel(gr_vars.behav_id((gr_vars.sleep == i-1 & gr_vars.insight == 0))))
end

sprintf('group N2 long: insight.yes: %d, no: %d', ...
numel(gr_vars.behav_id((gr_vars.sleep == 2 & gr_vars.n2_long == 1 & gr_vars.insight == 1))),...
numel(gr_vars.behav_id((gr_vars.sleep == 2 & gr_vars.n2_long == 1 & gr_vars.insight == 0))))

%% Fishers exact test 
gr1_1 = numel(gr_vars.behav_id((gr_vars.sleep == 2 & gr_vars.n2_long == 1 & gr_vars.insight == 1)));
gr1_2 = numel(gr_vars.behav_id((gr_vars.sleep == 2 & gr_vars.n2_long == 1 & gr_vars.insight == 0)));

gr2_1 = numel(gr_vars.behav_id((gr_vars.sleep == 1 & gr_vars.insight == 1)));
gr2_2 = numel(gr_vars.behav_id((gr_vars.sleep == 1 & gr_vars.insight == 0)));

tbl_ = [gr1_1, gr1_2; ...
    gr2_1, gr2_2];

[~, p, stats_] = fishertest(tbl_)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 3A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% please note that figure 3A and 3D were additionally plotted in R for consistency across the paper
% but numerically, plots here and in R are identical

contrast_ = 'sleep'; % 'insight' 
chOI = 'C4';

sleepsight_cluster_plot(dat, gr_vars, gitdir, contrast_, chOI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% models across channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract slope values per participant/channel to run stats model on each channel
restoredefaultpath % to show pvalues in fitglm of betas

gr_vars_table = table;
gr_vars_table.insight = gr_vars.insight';
gr_vars_table.sleep = num2str(gr_vars.sleep'); % make sleep a factor

for isub = 1:numel(dat)
    sl_(isub,:) = dat{isub}.powspctrm; % subj x channels
end

%% make it negative
sl_ = sl_ * -1;

%% specification of models

m01 = 'insight ~ 1';
m02 = 'insight ~ 1 + sleep';
m03 = 'insight ~ 1 + sleep + sl';
m04 = 'insight ~ 1 + sl';

%% lets go

chOI = {'F4', 'C4'}; % just the results for this channel will be printed

for ich = 1:size(dat{1}.powspctrm,1)

    gr_vars_table.sl = sl_(:,ich);

    m01_res = fitglm(gr_vars_table, m01, 'Distribution', 'binomial');
    m02_res = fitglm(gr_vars_table, m02, 'Distribution', 'binomial');
    m03_res = fitglm(gr_vars_table, m03, 'Distribution', 'binomial');
    m04_res = fitglm(gr_vars_table, m04, 'Distribution', 'binomial');

    
    % this prints stuff from the results text section (same order as it is reported)
    if sum(strcmpi(dat{1}.label{ich}, chOI)) == 1 % if one of the chOIs is the current channel, PRINT

        sprintf('------------- ------------- ------------- \n------------- results channel %s ------------- ', ...
            dat{1}.label{ich})
        sprintf('AIC. null: %2.1f AIC.sleep stage: %2.1f', dat{1}.label{ich}, ...
            m01_res.ModelCriterion.AIC, m02_res.ModelCriterion.AIC)
        sprintf('likelihood ratio null vs. sleep stage: %2.2f, p = %2.2f', m02_res.devianceTest.chi2Stat(2),...
            m02_res.devianceTest.pValue(2))

        % sleep stage model
        m02_res
        % full model AIC
        sprintf('AIC. sleep stage: %2.1f AIC.slope+sleep stage: %2.1f', m02_res.ModelCriterion.AIC, m03_res.ModelCriterion.AIC)

        % likelihood ratio test between models 
        lr = abs(2*(m03_res.LogLikelihood - m02_res.LogLikelihood));
        m03_res
        sprintf('likelihood ratio sleep stage vs. slope+sleep stage: %2.2f, p = %2.2f', lr,...
            1 - chi2cdf(lr, 1))

        % slope model AIC
        sprintf('AIC. slope+sleep stage: %2.1f AIC. slope: %2.1f', m03_res.ModelCriterion.AIC, m04_res.ModelCriterion.AIC)
        m04_res
        lr = abs(2*(m03_res.LogLikelihood - m04_res.LogLikelihood));
        sprintf('likelihood ratio slope+sleep stage vs. slope: %2.2f, p = %2.2f', lr,...
            1 - chi2cdf(lr, 2))
    end


    % MCFadden + AIC differences for topoplot (3B)

    % this calculates the adjusted McFadden pseudo Rsquare (penalises number of predictors)
    %m02_pseudo_r_2_adj(ich) = 1 - ((m02_res.LogLikelihood - (length(m02_res.Coefficients.Estimate) - 1)) / m01_res.LogLikelihood);
    m03_pseudo_r_2_adj(ich) = 1 - ((m03_res.LogLikelihood - (length(m03_res.Coefficients.Estimate) - 1)) / m01_res.LogLikelihood);
    m04_pseudo_r_2_adj(ich) = 1 - ((m04_res.LogLikelihood - (length(m04_res.Coefficients.Estimate) - 1)) / m01_res.LogLikelihood);
    

    % significant model fit (vs. null model)
    m03_pval(ich) = m03_res.devianceTest.pValue(2);
    m04_pval(ich) = m04_res.devianceTest.pValue(2);

    % difference in AIC scores
    AIC_diff03_02(ich) = m03_res.ModelCriterion.AIC - m02_res.ModelCriterion.AIC;
    AIC_diff04_03(ich) = m04_res.ModelCriterion.AIC - m03_res.ModelCriterion.AIC;


    


end


%% plot

addpath(genpath(gitdir)) % has to be added for color map (viridis)

varOI = AIC_diff04_03;
valOI = m04_pseudo_r_2_adj;
pvalOI = m04_pval;
label_idx = 2;

axislabel_ = {'R^2_{adjMcF} [full model]', 'AIC [full > baseline model]'; ...
    'R^2_{adjMcF} [slope model]', 'AIC [slope > full model]'};
    

idx_ = ismember(dat{1}.label, dat{1}.label(varOI < 0)) & ...
    (pvalOI <= 0.05)';
chOI_ = dat{1}.label(idx_); 
new_caxis = [0 0.083];
new_caxis_label = axislabel_{label_idx, 1};

col_map = viridis; 
tmp = struct;
tmp.freq = 1;
tmp.dimord = 'chan_freq';
tmp.label = dat{1}.label;
tmp.plotme  = varOI;
c_axis_range = [-3.2, 3.2];


figure; set(gcf, 'position', [0 0 1500 1500])
topo = [];
topo.layout = layout;
topo.parameter = 'plotme';
topo.gridscale = 360;
topo.marker = 'off';
topo.comment = 'no';
topo.style = 'straight';
cfg.shading = 'interp';

topo.highlight = 'on';
topo.highlightchannel = chOI_;
topo.highlightsymbol = '.';
topo.highlightsize = 100;
topo.highlightcolor = 'white';

ft_topoplotER(topo, tmp);
colormap(col_map)
set(gcf, 'Color', 'w');
hcb = colorbar;
hcb.Label.String = axislabel_{label_idx, 2};
hcb.FontSize = 40;
hcb.Label.FontSize = 60;
clim(c_axis_range);


% Highlight specific channels
hold on; % Keep the current plot
% getting corresponding values
idx_p = cellfun(@(x) find(strcmp(layout.label, x)), chOI_, 'UniformOutput', false);
idx_p = cell2mat(idx_p);
vals_ = valOI(idx_p);
new_col_amp = flipud(magma(256));
highlight_colors = interp1(linspace(new_caxis(1), new_caxis(2), size(new_col_amp, 1)), new_col_amp, vals_);

for ich = 1:numel(chOI_)
    chan_pos = layout.pos(idx_p(ich),:);
    scatter(chan_pos(1), chan_pos(2), 500, highlight_colors(ich, :), 'filled'); 
end

% Add a custom colorbar for the graded highlights
c = colorbar('westoutside'); % Add a colorbar for the graded highlights
colormap(c, new_col_amp); % Set the custom colormap
c.Limits = [-4, 4]; % Set the range for the graded highlights
c.TickLabels = round(linspace(new_caxis(1), new_caxis(2), 5),2); % Set the range for the graded highlights
c.Label.String = new_caxis_label; % Label the colorbar
c.Label.Interpreter = 'tex'; % Enable LaTeX interpreter
c.FontSize = 40;
c.Label.FontSize = 60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 3C & 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% please note that figure 3A and 3D were additionally plotted in R for consistency across the paper (see script ...) 
% but numerically, plots here and in R are identical
addpath(genpath(gitdir))

contrast_ = 'insight'; % 'insight' 
chOI = 'C4';

sleepsight_cluster_plot(dat, gr_vars, gitdir, contrast_, chOI)


%% save for plotting in R

chOI = 'C4';
chOI_idx = find(ismember(dat{1}.label, chOI));

slope_ = sl_(:,chOI_idx);

out_ = array2table([gr_vars.eeg_id',...
    gr_vars.sleep',...
    gr_vars.insight',...
    slope_]);
out_.Properties.VariableNames(1:4) = {'eeg_id', ...
    'sleep', ...
    'insight', ...
    'slope'};

savename = '1fslope_C4.csv';

if ~exist(fullfile(savedir, savefolder))
    mkdir(fullfile(savedir, savefolder))
end
writetable(out_, fullfile(savedir, savefolder, savename)) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% supplement table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath % to show pvalues in fitglm of betas

chOI = {'Fz', 'F3', 'F4', ...
    'Cz', 'C3', 'C4', ...
    'Pz', 'Oz'};


for ich = 1:numel(chOI)

    i = find(ismember(dat{1}.label, chOI{ich}));

    gr_vars_table.sl = sl_(:,i);

    m01 = 'insight ~ 1';
    m03 = 'insight ~ 1 + sleep + sl';
    m04 = 'insight ~ 1 + sl';

    m01_res = fitglm(gr_vars_table, m01, 'Distribution', 'binomial');
    m03_res = fitglm(gr_vars_table, m03, 'Distribution', 'binomial');
    m04_res = fitglm(gr_vars_table, m04, 'Distribution', 'binomial');

    % this calculates the adjusted McFadden pseudo Rsquare (penalises number of predictors)
    m03_pseudo_r_2_adj = 1 - ((m03_res.LogLikelihood - (length(m03_res.Coefficients.Estimate) - 1)) / m01_res.LogLikelihood);
    m04_pseudo_r_2_adj = 1 - ((m04_res.LogLikelihood - (length(m04_res.Coefficients.Estimate) - 1)) / m01_res.LogLikelihood);

    sprintf(['ch: %s\n\n' ...
        'Chi: %2.2f. p: %2.4f. Rsquared: %2.4f\n' ...
        'slope: %2.2f. p: %2.4f\n' ...
        'AIC: %2.2f\n\n' ...
        'Chi: %2.2f. p: %2.4f. Rsquared: %2.4f\n' ...
        'slope: %2.2f. p: %2.4f\n' ...
        'AIC: %2.2f\n'], ...
        dat{1}.label{i}, ...
        m03_res.devianceTest.chi2Stat(2), ...
        m03_res.devianceTest.pValue(2), ...
        m03_pseudo_r_2_adj, ...
        m03_res.Coefficients.Estimate(4), ...
        m03_res.Coefficients.pValue(4), ...
        m03_res.ModelCriterion.AIC, ...
        m04_res.devianceTest.chi2Stat(2), ...
        m04_res.devianceTest.pValue(2), ...
        m04_pseudo_r_2_adj, ...
        m04_res.Coefficients.Estimate(2), ...
        m04_res.Coefficients.pValue(2), ...
        m04_res.ModelCriterion.AIC)

end

