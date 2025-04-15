function sleepsight_cluster_plot(dat, gr_vars, gitdir, compOI, chOI)


load(fullfile(gitdir, 'mpib-60Ch-neighbours.mat')) %neighbours
load(fullfile(gitdir, 'mpib-60Ch-layout.mat')) %layout

%%
if nargin < 4
    compOI = 'insight'; % 'sleep'
    chOI = 'C4';
end

if nargin < 5
    chOI = 'C4';
end

%% SLEEP
cfg = [];
cfg.keepindividual = 'yes';

ga_wake = ft_freqgrandaverage(cfg, dat{gr_vars.sleep == 0});
ga_n1 = ft_freqgrandaverage(cfg, dat{gr_vars.sleep == 1});
ga_n2 = ft_freqgrandaverage(cfg, dat{gr_vars.sleep == 2});

%% INSIGHT
cfg = [];
cfg.keepindividual = 'yes';

ga_insight_yes = ft_freqgrandaverage(cfg, dat{gr_vars.insight == 1});
ga_insight_no = ft_freqgrandaverage(cfg, dat{gr_vars.insight == 0});

%%
ga = {};

if strcmpi(compOI, 'sleep')
    ga{1} = ga_wake;
    ga{2} = ga_n1;
    ga{3} = ga_n2;

elseif strcmpi(compOI, 'insight')
    ga{1} = ga_insight_yes;
    ga{2} = ga_insight_no;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stats. cluster based perm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg  = [];
%cfg.frequency = settings_.FOI;
cfg.frequency = 1;
cfg.channel = 'all';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.tail = 0;
cfg.alpha = .025; % two sided correction

if numel(ga) == 3
    cfg.statistic = 'ft_statfun_indepsamplesF';
    cfg.tail = 1; % For a F-statistic, it only make sense to calculate the right tail
    cfg.alpha = .05; % not two sided
end

cfg.method = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm = 'cluster'; % 'no', cluster;

cfg.clusteralpha = .05;
cfg.clustertail = cfg.tail;
cfg.correcttail = 'alpha'; % alpha prob no
cfg.neighbours = neighbours;
cfg.minnbchan = 2;
cfg.avgovertime = 'no';
cfg.avgoverchan = 'no';
cfg.avgoverfreq = 'no';
cfg.computecritval = 'yes';
cfg.numrandomization = 1000;
cfg.clusterstatistic = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'

nsub1 = size(ga{1}.powspctrm,1);
nsub2 = size(ga{2}.powspctrm,1);
nsub = nsub1 + nsub2;

if numel(ga) == 3
    nsub3 = size(ga{3}.powspctrm,1);
    nsub = nsub1 + nsub2 + nsub3;
    design = zeros(1,nsub);
    design(1:nsub1) = 1;
    design(nsub1+1:nsub1+nsub2) = 2;
    design(nsub1+nsub2+1:nsub1+nsub2+nsub3) = 3;
else
    design = zeros(1,nsub);
    design(1:nsub1) = 1;
    design(nsub1+1:nsub1+nsub2) = 2;
end

cfg.design = design;
cfg.ivar = 1;

% run stats
if numel(ga) == 3
    [Fieldtripstats] = ft_freqstatistics(cfg, ga{1}, ga{2}, ga{3});
else
    [Fieldtripstats] = ft_freqstatistics(cfg, ga{1}, ga{2});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot plot plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%% TOPOS
%%%%%%%%%%%%%%%%%

%% topoplot just for INSIGHT

if strcmpi(compOI, 'insight')

    % to plot data
    dat_insight = mean(ga_insight_yes.powspctrm,1);
    dat_noinsight = mean(ga_insight_no.powspctrm,1);

    dat_diff = dat_insight*-1 - dat_noinsight*-1;

    prob_thres = .05; % threshold that has to be reached for plotting

    col_amp = viridis;

    for i = 1:numel({Fieldtripstats.posclusters.prob})

        if Fieldtripstats.posclusters(i).prob < prob_thres
            b = Fieldtripstats.stat .* (Fieldtripstats.posclusterslabelmat==i);
            tsums = squeeze(nansum(b,2));
            tsums = repmat(tsums,[1 size(Fieldtripstats.stat,2)]);

            sigmat = Fieldtripstats.posclusterslabelmat==i;
            sigfreq = round(Fieldtripstats.freq(any(sigmat)));

            tmp = Fieldtripstats;
            %tmp.plotme  = Fieldtripstats.stat;
            tmp.plotme  = dat_diff;
            %tmp.plotme  = tsums;

            figure; set(gcf, 'position', [0 0 1500 1500])
            topo = [];
            topo.layout = layout;
            topo.parameter = 'plotme';
            topo.gridscale = 360;
            topo.marker = 'off';
            topo.comment = 'no';
            topo.style = 'straight';
            %topo.shading = 'interp';

            topo.highlight = 'on';
            topo.highlightchannel = Fieldtripstats.label(Fieldtripstats.posclusterslabelmat == 1);
            topo.highlightsymbol = '.';
            topo.highlightsize = 40;
            topo.highlightcolor = 'white';

            ft_topoplotER(topo, tmp);
            colormap(col_amp)
            caxis([-0.3,0])
            set(gcf, 'Color', 'w');

            hcb = colorbar;
            hcb.FontSize = 40;
            hcb.Ticks = -0.3:0.1:0;
            hcb.Label.String = {'spectral slope', '[Insight > No Insight]'};
            hcb.Label.FontSize = 60;

        end
    end

end
%% 

chOI_idx = find(ismember(dat{1}.label, chOI));

for isub = 1:nsub

    sl_(isub) = dat{isub}.powspctrm(chOI_idx) *-1;

end

%%
font_size = 20;

if strcmpi(compOI, 'insight') == 1

    font_size = 20;
    pow1 = sl_(gr_vars.insight == 1);
    pow2 = sl_(gr_vars.insight == 0);

    g1 = repmat({'insight'},numel(pow1),1);
    g2 = repmat({'no'},numel(pow2),1);
    g = [g1; g2];
    figure()
    boxplot([pow1, pow2], g, 'orientation', 'horizontal')
    hold on
    plot([mean(pow1), mean(pow2)], [1:2], '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 12)

    ax = gca;
    ax.FontSize = font_size;
    ylabel({'1/f slope'}, FontSize = font_size)

    sprintf('mean 1/f slope. insight = %2.2f. %2.2f.\n no insight = %2.2f. %2.2f', mean(pow1), std(pow1)/sqrt(numel(pow1)), ...
        mean(pow2), std(pow2)/sqrt(numel(pow2)))

    [~, p, ~, stat_] = ttest2(pow1, pow2, 'Vartype', 'unequal');
    sprintf('in vs. no in: t = %2.2f. p = %2.2f. d = %2.2f', stat_.tstat, p, computeCohen_d(pow2, pow1))

elseif strcmpi(compOI, 'sleep') == 1

    pow1 = sl_(gr_vars.sleep == 0);
    pow2 = sl_(gr_vars.sleep == 1);
    pow3 = sl_(gr_vars.sleep == 2);

    g1 = repmat({'wake'},numel(pow1),1);
    g2 = repmat({'n1'},numel(pow2),1);
    g3 = repmat({'n2'},numel(pow3),1);
    g = [g1; g2; g3];
    figure()
    boxplot([pow1, pow2, pow3], g, 'orientation', 'horizontal')
    hold on
    plot([mean(pow1), mean(pow2), mean(pow3)], [1:3], '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 12)

    ax = gca;
    ax.FontSize = font_size;
    xlabel({'1/f slope'}, FontSize = font_size)

    sprintf('mean and sem 1/f slope. \nw = %2.2f. %2.2f.\nn1 = %2.2f. %2.2f.\nn2 = %2.2f. %2.2f.', mean(pow1), std(pow1)/sqrt(numel(pow1)), ...
        mean(pow2), std(pow2)/sqrt(numel(pow2)), ...
        mean(pow3), std(pow3)/sqrt(numel(pow3)))


    %% t-tests for channel x
    [~, p, ~, stat_] = ttest2(pow1, pow2, 'Vartype', 'unequal');
    sprintf('w vs. n1: t = %2.2f. p = %2.2f. d = %2.2f', stat_.tstat, p, computeCohen_d(pow1, pow2))
    
    [~, p, ~, stat_] = ttest2(pow2, pow3, 'Vartype', 'unequal');
    sprintf('n1 vs. n2: t = %2.2f. p = %2.2f. d = %2.2f', stat_.tstat, p, computeCohen_d(pow2, pow3))

end



