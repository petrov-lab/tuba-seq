import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def percentile_plot(data, ax, order, 
    baseline=True, #percentiles=None, 
    hue_map=None, alpha=0.05, inert_darkness=0.25, sgRNA_spacing=0.1, saturation_of_lightest_percentile=1, palette_constructor=sns.light_palette,
    add_legend=True, legend_text_size='x-small',
    dot_kwargs=dict(marker='.', linewidths=1, edgecolors='k', s=400, zorder=3),
    ebar_kwargs=dict(ecolor='#404040', elinewidth=0.5, ls='none', capsize=400, zorder=2),
    xtick_kwargs=dict(rotation=90, style='italic', size=20, ha='center'),
    baseline_kwargs=dict(color='k', ls='dotted', lw=1)):

    if hue_map is None:     # Generate unique hues for each target, if not provided.
        hue_map = sns.husl_palette(n_colors=len(order), s=1, l=inert_darkness) 
    
    for i, target in enumerate(order):        # Plot percentile curves for each sgID 1-by-1
        df = data.query('target == @target')
        X = i + np.linspace(0, 1-sgRNA_spacing, num=len(df), endpoint=False) 
        Y = df['true']
        
        err = np.vstack((df['true'] - df['low'], df['high'] - df['true']) )
        if err.any():
            ax.errorbar(X, Y, yerr=err, **ebar_kwargs)
        
        n_colors = len(Y) + saturation_of_lightest_percentile
        inert_colors = palette_constructor(3*(inert_darkness,), n_colors=n_colors)[saturation_of_lightest_percentile:]
        hue_colors = palette_constructor(hue_map[target], n_colors=n_colors)[saturation_of_lightest_percentile:]
        colors = [hue if pval < alpha else inert for hue, inert, pval in zip(hue_colors, inert_colors, df['P-value'])]
        ax.scatter(X, Y, c=colors, label=target, **dot_kwargs)
    
    if baseline:
        ax.axhline(1, **baseline_kwargs)
    
        #Y-axis specifications
    ax.set_yscale('log', basey=2)
    from matplotlib import ticker
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.set_ylabel("Relative Size (sg$TS$/sg$Inerts$)") 
    
        #X-axis specifications
    ax.xaxis.set_ticks(np.arange(len(order))+(1-sgRNA_spacing)/2)
    ax.xaxis.set_ticklabels(order, **xtick_kwargs)
    if add_legend:
        X = np.linspace(sgRNA_spacing, 1, num=len(df)+1)
        y = ax.get_ylim()[1]*0.99   # hack to get top of plot 
        percentiles = df.index.get_level_values('Percentile').values.tolist()
        ax.scatter(X[:-1], len(percentiles)*[y], c=inert_colors, label='Legend', **dot_kwargs)
        for x, p in zip(X, percentiles):
            ax.text(x, y, '{:g}  '.format(p), ha='center', va='top', rotation=90, size=legend_text_size)
        ax.text(X[-1], y,'Percentile', va='center', ha='left', size=legend_text_size)
    return ax

def poisson_CI(counts, alpha=0.05, min_high_CI=1.0):
    from scipy.stats.distributions import poisson
    CI = pd.DataFrame(list(poisson.interval(0.95, [4, 3, 1])), index=['low', 'high']).T
    CI = CI.fillna(0)
    CI['high'] = CI['high'].clip(min_high_CI)
    return CI

def errorbar_line_histogram(X, bins=14, error_style='shadow', error_func=poisson_CI, shadow_alpha=0.28, points=False, y_floor=1e-16, 
                            normed=False, trim_flanking_zeros=True, xscale='log', yscale='log', **kargs):
    N = len(X)
    counts, bin_edges = np.histogram(np.array(X), bins=bins)
    dx = np.diff(bin_edges)
    density = counts/(N*dx)
    CI = error_func(counts)
    
    ax = kargs.pop('ax', plt.gca())
    if xscale == 'log':
        log_edges = np.log(bin_edges)
        X = np.exp(log_edges + np.diff(log_edges)/2)
    else:
        X = bin_edges[:-1] + dx/2
   
    if normed:
        Y = density 
        CI /= N*dx
    else:
        Y = counts

    CI['low'] = CI['low'].clip(y_floor)
    
    line, = ax.plot(X, Y, **kargs)
    color = line.get_color()
    
    if error_style == 'shadow':
        ax.fill_between(X, CI['low'], CI['high'], edgecolor='none', facecolor=color, alpha=shadow_alpha)
    elif error_style == 'bar':
        ax.errorbar(X, Y, yerr=(CI[['low', 'high']] - Y).abs().values.T, color=color, ecolor=color, elinewidth=2)
    else:
        print('A valid error_style was not given.')
    ax.set(xscale=xscale, yscale=yscale)
    return counts, bins, line

def LN_PL_best_fit_pdf(S, PL_summary, ax, label, color, 
                        poor_fit_alpha=0.2, xmin=3, decades=4, 
                        double_log_bins=True, double_log_offset=1.5, bins=14, resolution=5000, 
                        alpha=0.05, ylim=[1e-12, 1e-3], **kargs):
    
    X = np.logspace(xmin, xmin+decades, resolution)
    xlim = X[[0, -1]]
    if double_log_bins:
        final_bins = pow(10, np.logspace(*np.log10(double_log_offset+np.array([0, decades])), num=bins))
        final_bins *= xlim[0]/final_bins[0]
    else:
        final_bins = np.logspace(xmin, xmin+decades, bins)

    kwargs = dict(color=color, label=label) 
    counts, bins, line = errorbar_line_histogram(S, ax=ax, bins=final_bins, normed=True, xscale='log', yscale='log', ls='o', lw=3, **kwargs)
    
    #mean = S.mean()
    #var = S.var()
    #mu = np.log(mean**2/np.sqrt(var+mean**2))
    #s = np.sqrt(np.log(1 + var/mean**2))
    #N = X*s*np.sqrt(2*np.pi)
    #LN_pdf = np.exp(-(np.log(X)-mu)**2/(2*s**2))/N

    from scipy.stats.distributions import lognorm
    LN_pdf = lognorm.pdf(X, lognorm.fit(S, floc=0))
   
    kwargs['lw'] = 2
    kwargs['ls'] = '-'
    
    if PL_summary['Best Fit'] == 'LN':
        ax.plot(X, LN_pdf, **kwargs)
    else:
        x_min = PL_summary['x_min']
        PL_range = X >= x_min
        X_PL = X[PL_range]
        
        ax.plot(X[-PL_range], LN_pdf[-PL_range], **kwargs)
        ax.plot(X_PL, LN_pdf[PL_range], alpha=poor_fit_alpha, **kwargs)
            # Try scipy's thing 
        a = PL_summary['alpha']
        Y_PL = a*x_min**(a - 1)*X_PL**-a
        PL_fraction = (S > x_min).mean()
        Y_PL *= PL_fraction

        ax.plot(X_PL, Y_PL, **kwargs)

    ax.text(xlim[1], ylim[1], r"""$P = {P-value:.4f}$
$\alpha = {alpha:.2f}$""".format(**PL_summary), ha='right', va='top')
    ax.set( ylabel='Probability Density', xlabel='Cells', title=label, xlim=xlim, ylim=ylim )
    return ax

def fancy_percentage_formatter(x, sig_figs=2):
    if x == 0.:
        return '0%'
    rx = round(x, -int(np.floor(np.log10(abs(x))) - (sig_figs - 1)))
    return '{:g}%'.format(rx*100)

def jitter_plot(S, order, colors=None, ax=plt.gca(),
                annotate_mean=True, tumor_numbers=True, decade_percentages=False, ymax=None,
                jitter=0.4, scale=5e-4, text_size='large', text_lw=3, mean_bar_width=0.9, 
                xtick_kwargs=dict(rotation=90, style='italic', size=20, ha='center'),
                mean_bar_kwargs=dict(color='k', lw=3.5, zorder=3)):
    
    import matplotlib.patheffects as path_effects 
    X = np.arange(len(order))
    xlim = -0.5, X[-1] + 0.5
    
    if colors is None:
        colors = dict(zip(order, sns.color_palette(n_colors=len(order))))
    
    if ymax is None:
        ymax = pow(10, np.ceil(np.log10(S.max()*(1.5 if tumor_numbers else 1))))

    gb = S.groupby(level='target')
    
    if annotate_mean:
        from tuba_seq.tools import LN_mean
        ax.hlines(gb.apply(LN_mean).loc[order], X-mean_bar_width/2, X+mean_bar_width/2, **mean_bar_kwargs)
    
    for X_i, rna in zip(X, order):
        Y = gb.get_group(rna).values
        x = X_i + 2*jitter*np.random.random(len(Y)) - jitter
        ax.scatter(x, Y, s=Y*scale, color=colors[rna], label=rna, zorder=10)
    
    if tumor_numbers:
        N = gb.count()[order].values
        for x, n in enumerate(N):
            ax.text(x, ymax, '$N=$\n${:,}$'.format(n), ha='center', va='top')
    if decade_percentages:
        decade_mins = pow(10, np.arange(np.floor(np.log10(S.min())), np.log10(S.max())))
        ax.hlines(decade_mins, *xlim, color='0.25', lw=1, linestyles='dashed')
        def decade_fraction(S, decade_min):
            return ((S >= decade_min)*(S < decade_min*10)).mean()
        
        for decade_min in decade_mins:
            fractions = gb.agg(decade_fraction, decade_min)
            if fractions.sum() == 0.:
                continue
            y = decade_min*np.sqrt(10)
            for x, frac in enumerate(fractions.loc[order].values):
                text = ax.text(x, y, fancy_percentage_formatter(frac), size=text_size, ha='center', va='center', color='w', weight='bold', zorder=15)
                text.set_path_effects([path_effects.Stroke(linewidth=text_lw, foreground='black', capstyle='round', joinstyle='round'), path_effects.Normal()])       
    
    ax.set(xlim=xlim, ylabel='Cells (Absolute no.)', yscale='log')
    ax.set_ylim(ymax=ymax)
    ax.xaxis.set_ticks(X)
    ax.xaxis.set_ticklabels(order, **xtick_kwargs)
    return ax

def de_step(X, Y, xlog=True):
    if xlog:
        X = np.log(X)
    Xmid = (X[::2] + X[1::2])/2
    return np.vstack((np.exp(Xmid) if xlog else Xmid, Y[::2])).T

def __errorbar_step_histogram(X, histtype='line', error='shadow', error_func=poisson_CI, error_transparency=0.28, points=False, y_floor=1e-16, bels=None, normed=False, trim_flanking_zeros=True, **kargs):
    ### BROKEN####
    assert histtype in ['line', 'step'], "Error Shading doesn't make sense in a barplot, histtype must be 'line' or 'step'."
    N = len(X)
    ax = kargs.pop('ax', plt.gca())
    density, bins, patches = ax.hist(np.array(X), histtype='step', normed=normed, **kargs)
    assert type(density) is np.ndarray, "errorbar_histogram only works with one distribution at a time"
        # Convert densities back to counts, if histogram is normalized
    dx = np.diff(bins)
    counts = np.round(density*N*dx if normed else density).astype(np.int64)
        # Get mid-point of each step--to use as a point  
    first_patch = patches[0]
    color = kargs.get('color', first_patch.get_facecolor())
    
    if trim_flanking_zeros:
        middle = slice(len(counts) - len(np.trim_zeros(counts, 'f')), len(np.trim_zeros(counts, 'b')))
        trimed = first_patch.get_xy()[1:-1]
    
    trimed[np.repeat(density==0, 2), 1] = y_floor
    if histtype == 'line':
        trimed = de_step(*trimed.T)
    X, Y = trimed.T
    if histtype == 'line':
        trimed = trimed[middle]

    first_patch.set_xy(trimed)
    if points:
        ax.plot(X[middle], Y[middle], '.', color=color, lw=3)
    if bels is not None:
        ymax = int(np.ceil(np.log10(Y).max()))
        ax.set_ylim(10**(ymax-bels), 10**ymax)
        
    CI = error_func(counts)
    # DE-INDENT
    if normed:
        CI /= N*dx
    CI['low'] = CI['low'].clip(y_floor)

    #if histtype == 'step':
    #    CI = np.repeat(low_CI, 2)
    #    high_CI = np.repeat(high_CI, 2)
    #if error == 'shadow':
    #    ax.fill_between(X, low_CI, high_CI, edgecolor='none', facecolor=color, alpha=transparency)
    #elif error == 'bar':
    #    ax.errorbar(X, Y, yerr=[Y - low_CI, high_CI - Y], color=color, ecolor=color, elinewidth=2)
    #else:
    #    print('A valid error choice was not given')
    return counts, bins, patches

class NullObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        pass

def text_color_legend(ax, visible_handles=False, legend_prop={'weight':'semibold'}, bbox_to_anchor=(1, 1), **kargs):
    """text_color_legend() -> eliminates legend key and simply colors labels with the color of the lines."""
    handles, labels = ax.get_legend_handles_labels()
    handles = [handle[0] if type(handles) == list else handle for handle in handles]
    if not visible_handles: 
        kargs['handler_map'] = {handle:NullObjectHandler() for handle in handles}
    L = ax.legend(handles, labels, prop=legend_prop, borderaxespad=0, bbox_to_anchor=bbox_to_anchor, **kargs)
    for handle, text in zip(handles, L.get_texts()):
        text.set_color(handle.get_facecolor() if handle.get_fill() else handle.get_edgecolor())
    return L



