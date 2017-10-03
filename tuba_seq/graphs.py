import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
import pandas as pd

def percentile_plot(data, ax, order, 
    baseline=True, percentiles=None, 
    hue_map=None, alpha=0.05, inert_darkness=0.25, sgRNA_spacing=0.1, saturation_of_lightest_percentile=1, palette_constructor=sns.light_palette,
    add_legend=True, legend_text_size='x-small',
    dot_kwargs=dict(marker='.', linewidths=1, edgecolors='k', s=350, zorder=3),
    ebar_kwargs=dict(ecolor='#404040', elinewidth=0.5, ls='none', capsize=400, zorder=2),
    xtick_kwargs=dict(rotation=90, style='italic', size=20, ha='center'),
    baseline_kwargs=dict(color='k', ls='dotted', lw=1)):

    if hue_map is None:     # Generate unique hues for each sgID, if not provided.
        hue_map = sns.husl_palette(n_colors=len(order), s=1, l=inert_darkness) 
    
    for i, sgID in enumerate(order):        # Plot percentile curves for each sgID 1-by-1
        df = data.query('sgID == @sgID')
        X = i + np.linspace(0, 1-sgRNA_spacing, num=len(df), endpoint=False) 
        Y = df['true']
        CI = df.loc[['low', 'high']] # extract the Confidence Interval manually. 
        err = (CI.T - Y.values).abs().values 
        
        if err.any():
            ax.errorbar(X, Y, yerr=err, **ebar_kwargs)
        
        n_colors = len(Y) + saturation_of_lightest_percentile
        inert_colors = palette_constructor(3*(inert_darkness,), n_colors=n_colors)[saturation_of_lightest_percentile:]
        hue_colors = palette_constructor(hue_map[sgID], n_colors=n_colors)[saturation_of_lightest_percentile:]
        colors = [hue if pval < alpha else inert for hue, inert, pval in zip(hue_colors, inert_colors, df['P-value'])]
        ax.scatter(X, Y, c=colors, label=sgID, **dot_kwargs)
    
    if baseline:
        ax.axhline(1, **baseline_kwargs)
    
        #Y-axis specifications
    ax.set_yscale('log', basey=2)
    from matplotlib.ticker import ScalarFormatter
    ax.yaxis.set_major_formatter(ScalarFormatter())
        #X-axis specifications
    ax.xaxis.set_ticks(np.arange(len(order))+(1-sgRNA_spacing)/2)
    ax.xaxis.set_ticklabels(order, **xtick_kwargs)
    
    if add_legend:
        X = np.linspace(sgRNA_spacing, 1, num=len(df), endpoint=False)
        y = ax.get_ylim()[1]*0.99   # hack to get top of plot 
        percentiles = df.index.get_level_values('Percentile').values.tolist()
        ax.scatter(X, len(percentiles)*[y], c=inert_colors, label='Legend', **dot_kwargs)
        for x, p in zip(X, percentiles):
            ax.text(x, y, p[:-1] + '  ', ha='center', va='top', rotation=90, size=legend_text_size)
        ax.text(X[-1], y,'  '+'Percentile', va='center', ha='left', size=legend_text_size)
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

def fancy_percentage_formatter(x):
    decimal_places = -np.floor(np.log10(x*100))
    return '{:g}%'.format(round(x*100, min(0, decimal_places)))

def jitter_plot(S, ax, order, colors,
                annotate_mean=True, tumor_numbers=True, decade_percentages=False,
                jitter=0.4, scale=5e-4, text_size='large', text_lw=3, mean_bar_width=0.9, 
                xtick_kwargs=dict(rotation=90, style='italic', size=20, ha='center'),
                mean_bar_kwargs=dict(color='k', lw=3.5, zorder=3)):
    
    import matplotlib.patheffects as path_effects 
    X = np.arange(len(order))
    xlim = -0.5, X[-1] + 0.5
    
    gb = S.groupby(level='sgID')
    
    if annotate_mean:
        from tools import LN_mean
        ax.hlines(gb.apply(LN_mean).loc[order], X-mean_bar_width/2, X+mean_bar_width/2, **mean_bar_kwargs)
    
    for i, rna in enumerate(order):
        Y = gb.get_group(rna).values
        X = i + 2*jitter*np.random.random(len(Y)) - jitter
        ax.scatter(X, Y, s=Y*scale, color=colors[rna], label=rna, zorder=10)
    
    if tumor_numbers:
        y = S.max()*1.05
        N = gb.count()[order].values
        for x, n in enumerate(N):
            ax.text(x, y, '$N=$\n${:,}$'.format(n), ha='center', va='bottom')
    if decade_percentages:
        sns.set_style('whitegrid')
        decade_mins = ax.get_yticks()[:-1]
        def decade_fraction(S, decade_min):
            return ((S >= decade_min)*(S < decade_min*10)).mean()
        
        for decade_min in decade_mins:
            fractions = gb.agg(decade_fraction, decade_min)
            if fractions.sum() == 0.:
                continue
            y = decade_min + np.sqrt(10)
            for x, frac in fractions.loc[order].values:
                text = ax.text(x, y, fancy_percentage_formatter(frac), size=text_size, ha='center', va='center', color='w', weight='bold', zorder=15)
                text.set_path_effects([path_effects.Stroke(linewidth=text_lw, foreground='black', capstyle='round', joinstyle='round'), path_effects.Normal()])       

    ax.set(xlim=xlim, ylabel='Cells (Absolute no.)', yscale='log')
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

def barcode_diversity_report(tumor_numbers, degeneracies=None, plot=True):
    sns.set_style('whitegrid')
        # Frequency profile
    barcodes = tumor_numbers.index.get_level_values('barcode')
    random_barcode_length = len(barcodes.values[0])
    PWM = pd.concat({i:barcodes.str.get(i).value_counts()/len(barcodes) for i in range(random_barcode_length)})
    PWM.index.names = ['Location', 'Nucleotide']
    PWM.name = 'Frequency'
    if plot:
        f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6, 12))
        base = np.zeros(random_barcode_length)
        for nuc in 'ACGT':
            X = PWM.loc[:, nuc]
            ax1.bar(X.index.values, X, bottom=base, label=nuc)
            base += X.values
        yticks = [0, 0.25, 0.5, 0.75, 1]
        ax1.set(ylabel='Frequency', xlabel='Position', ylim=[0, 1], yticks=yticks, yticklabels=list(map('{:.0%}'.format, yticks)))
        text_color_legend(ax1, title='Nucleotide:') 
        sns.despine(left=True, right=True, top=False)
    if degeneracies == None:
        print("Will infer the number of random bases from frequency profile of barcodes...")
        random_bases = (PWM > 0.1).groupby(level='Location').sum()
        base_degeneracies = random_bases.value_counts()
        base_degeneracies.index = np.int64(base_degeneracies.index.values)
        print("Found:")
        for k, v in base_degeneracies.items():
            if k != 1:
                print("{:} {:}-fold degenerate positions".format(v, k))
        
        multiplicity = base_degeneracies.index.values**base_degeneracies.values
        degeneracies = np.multiply.reduce(multiplicity)
        print("Thus, there are {:,} potential random barcodes".format(degeneracies))

    def sgID_statistics(S):
        M = S.unstack(level='Mouse').notnull()
        observed = len(M)
        barcodes_per_mouse = M.sum()
        N = len(barcodes_per_mouse)
        mu_geq1 = M.sum(axis=1).mean()
        from scipy.optimize import brentq
        def P_0(p):
            return (1 - p)**N
        p = brentq(lambda p: mu_geq1 - p*N/P_0(p), 0, mu_geq1/N, full_output=False)
        unobserved = observed*P_0(p)/(1 - P_0(p))

        def collision_probability(N):
            from scipy.stats.distributions import poisson
            P_0, P_1 = poisson.pmf([0, 1], N/(observed+unobserved))
            return (1 - P_0 - P_1)/(1 - P_0)
        
        P_collision_by_mouse = barcodes_per_mouse.apply(collision_probability)
        return {'Estimated Unobserved Barcodes' : unobserved,
                'Observed Barcodes' : observed,
                'Probability of Barcode Collision' : P_collision_by_mouse.dot(barcodes_per_mouse)/barcodes_per_mouse.sum()}

    sgID_info = tumor_numbers.groupby(level='target').apply(sgID_statistics).unstack()
    if plot:
        X = np.arange(len(sgID_info))
        bars = ax2.bar(X, sgID_info['Observed Barcodes'], label='Observed')
        ax2.bar(X, sgID_info['Estimated Unobserved Barcodes'], label='Estimated Unobserved',
                bottom=sgID_info['Observed Barcodes'], color=bars[0].get_facecolor(), alpha=0.5)
        ax2.set_ylabel('Number of Barcodes')
        text_color_legend(ax2, title='Barcodes:', bbox_to_anchor=(0.9, 1)) 
        
        ax3.bar(X, sgID_info['Probability of Barcode Collision'])
        max_ytick = int(np.ceil(sgID_info['Probability of Barcode Collision'].max()*100))
        yticks = np.arange(max_ytick+1)/100
        ax3.set(ylabel='Probability of Barcode Collision', yticks=yticks, ylim=[0, max_ytick/100], yticklabels=list(map('{:.0%}'.format, yticks)))
        for ax in (ax2, ax3):
            ax.set(xticks=X, xlabel='Target')
            sns.despine(left=True, right=True, top=False, ax=ax)
            ax.xaxis.set_ticklabels(sgID_info.index.values, rotation=90, style='italic')
            ax.xaxis.grid(False)

        plt.savefig('barcode_diversity_report.pdf', format='pdf', bbox_inches='tight')
    return sgID_info