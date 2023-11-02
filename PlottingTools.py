# Author: Alex Baras, MD, PhD (https://github.com/alexbaras)
# NCATS Maintainer: Dante J Smith, PhD (https://github.com/djsmith17)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import altair as alt
alt.data_transformers.disable_max_rows()
from scipy import ndimage as ndi


def plot_2d_density(X, Y=None, bins=200, n_pad=40, w=None, ax=None, gaussian_sigma=0.5, cmap=plt.get_cmap('viridis'), vlim=np.array([0.001, 0.98]), circle_type='bg', box_off=True, return_matrix=False):
    '''plot_2d_density(X, Y, bins, n_pad, w, ax, gaussian_sigma, cmap, vlim, circle_type, box_off, return_matrix)
    is a method for drawing 2D histograms figures. In this particular instance, we are plotting the outputs of the UMAP.
    
    Parameters:
    X: X-values
    Y: Y-values
    bins: Bins. Vector tickMin and tickMax for X and Y direction.
    n_pad: number of elements to pad around a bin to make the spectrogram a little prettier. (I think)
    w: weights used to augment the density maps based on secondary phenotypes
    ax: Axes handle for the current figure
    gaussian_sigma: Sigma value for a gaussian filter with weight values
    cmap: color-map to use for the 2D Histogram
    vlim: Value range that the colormap should cover
    circle_type: Type of colormesh to return. ('bg', 'arch'). Will need to investigate further
    box_off: Flag to turn box around the figure off
    return_matrix: Flag for returning matrix. Default = 0 to not return density matrix, but intstead draws the figure.

    Returns:
    
    '''
    if Y is not None:
        if w is not None:
            b, _, _ = np.histogram2d(X, Y, bins=bins)
            b = ndi.gaussian_filter(b.T, sigma=gaussian_sigma)

            s, _, _ = np.histogram2d(X, Y, bins=bins, weights=w)
            s = ndi.gaussian_filter(s.T, sigma=gaussian_sigma)

            d = np.zeros_like(b)
            # d[b > 0] = s[b > 0] / b[b > 0]
            d = s
            d = ndi.gaussian_filter(d, sigma=gaussian_sigma)
        else:
            d, _, _ = np.histogram2d(X, Y, bins=bins)
            d /= np.sum(d)
            d = ndi.gaussian_filter(d.T, sigma=gaussian_sigma)
    else:
        d = X

    if return_matrix:
        return d
    else:
        if np.isscalar(vlim):
            vlim = np.array([0, np.quantile(d[d > 0].flatten(), vlim)])
        else:
            if np.all((vlim < 1) & (vlim > 0)):
                vlim = np.quantile(d[d > 0].flatten(), vlim)

        if ax is None:
            _, ax = plt.subplots()

        if np.isscalar(bins):
            n_bins = bins
        else:
            n_bins = len(bins[0]) - 1

        if circle_type == 'bg':
            c = np.meshgrid(np.arange(2 * n_pad + n_bins), np.arange(2 * n_pad + n_bins))
            c = np.sqrt(((c[0] - ((2 * n_pad + n_bins) / 2)) ** 2) + ((c[1] - ((2 * n_pad + n_bins) / 2)) ** 2)) < (0.95 * ((2 * n_pad + n_bins) / 2))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]) + c, vmin=1, vmax=1 + vlim[1], cmap=cmap, shading='gouraud', alpha=1)
            # ax.pcolormesh(np.log10(np.pad(d, [n_pad, n_pad]) + c + 1), vmin=np.log10(2), vmax=np.log10(2 + vlim[1]), cmap=cmap, shading='gouraud', alpha=1)
        elif circle_type == 'arch':
            c = (n_bins / 2)
            ax.add_artist(plt.Circle((c + n_pad, c + n_pad), 0.95 * (c + n_pad), color='black', fill=False))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=np.quantile(d[d < 0].flatten(), 0.03), 
                                                     vmax=np.quantile(d[d > 0].flatten(), 0.97), cmap=cmap, shading='gouraud', alpha=1)
        else:
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=0, vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)

        if box_off is True:
            [ax.spines[sp].set_visible(False) for sp in ax.spines]
            ax.set(xticks=[], yticks=[])


def plt_cmap(ax, cmap, extend, width, ylabel):
    '''plt_cmap(ax, cmap, extend, width, ylabel) draws a colorbar for the current colormap at the correct
    axes location, and with the correct label.

    Parameters:
    ax: Matplotlib axes handle
    cmap: Matplotlib colormap
    extend: {'neither', 'both', 'min', 'max'}  
            Make pointed end(s) for out-of-range values (unless 'neither'). 
            These are set for a given colormap using the colormap set_under and set_over methods.
    width: Width of the colorbar in Figure coordinates. '0.01' suggested value
    ylabel: String of the colomap label

    Returns:

    '''
    cb = mpl.colorbar.ColorbarBase(ax=ax, cmap=cmap, extend=extend)
    cb.set_ticks([])
    pos = ax.get_position().bounds
    ax.set_position([pos[0], pos[1], width, pos[3]])
    ax.set(ylabel=ylabel)


def plot_spatial_elem(ax, elems, title, color): 
    '''Generates a scatter plot of the cell positions (X/Y) from the sample collected
    '''
    ax.set_title(title)
    ax.set_xlabel('Centroid X')
    ax.set_ylabel('Centroid Y')
    spatialMax = [elems['X0'].unique()*2, elems['Y0'].unique()*2]
    ax.set_xlim([0, spatialMax[0]])
    ax.set_ylim([0, spatialMax[0]])
    ax.set_frame_on(False)
    plt.scatter(elems['CentroidX'], elems['CentroidY'], s=2, c=color)


def plot_spatial_interactive(elems, title, feature):
    '''Generates a scatter plot of the cell positions (X/Y) from the sample collected
    '''
    selection = alt.selection_multi(fields=[feature], bind='legend')
    chart = alt.Chart(elems).mark_circle().encode(
            x='CentroidX',
            y='CentroidY',
            color= alt.Color(feature, legend=alt.Legend(
                                                            orient='bottom',
                                                            columns = 4)),
            opacity=alt.condition(selection, alt.value(1), alt.value(0.2)),
            tooltip=feature
            ).properties(width=750,height=750
            ).interactive().add_selection(selection)
    
    return chart


def plot_neighborhood_profile(ax, cell_label, dist_bin, cell_density, phenoSet, maxDens=0.1, legF=0):
    '''This function generates the line plots of the phenotype density 
    at different distances from a given cell
    '''

    SlBgC  = np.array([14, 17, 23])/256    # Streamlit Background Color
    SlTC   = np.array([250, 250, 250])/256 # Streamlit Text Color
    Sl2BgC = np.array([38, 39, 48])/256    # Streamlit Secondary Background Color

    plotaxes = []
    plotLabels = []
    for i , key in enumerate(phenoSet):
        plotax, = ax.plot(dist_bin, cell_density[:, i], color=phenoSet[key])
        ax.fill_between(dist_bin, cell_density[:, i], color=phenoSet[key], alpha=.35)

        ax.set_xticks(dist_bin)
        # ax.set_xticklabels(['0-25', '26-50', '51-100', '101-150', '151-200'])
        ax.set_xlim([25, 200])
        ax.set_ylim([0, maxDens])
        ax.set_title(f'Cell #{cell_label}', fontsize = 16, color = SlTC)
        ax.set_xlabel('Spatial Bound (\u03BCm)', fontsize = 14, color = SlTC)
        ax.set_ylabel('Cell Density', fontsize = 14, color = SlTC)

        ax.set_frame_on(False)
        ax.spines['left'].set_color(SlTC)
        ax.spines['bottom'].set_color(SlTC)
        ax.tick_params(axis='x', colors=SlTC, which='both')
        ax.tick_params(axis='y', colors=SlTC, which='both')
        plotaxes.append(plotax)
        plotLabels.append(key)

    if legF:
        ax.legend(plotaxes, plotLabels,
                bbox_to_anchor=(-0.05, -0.1), 
                loc='upper left', 
                borderaxespad=0, 
                ncols = 4,
                facecolor = Sl2BgC,
                edgecolor = Sl2BgC,
                labelcolor = SlTC)


def plot_neighborhood_profile_propor(ax, cell_label, dist_bin, cell_propor, phenoSet, colors, legF=0):
    '''This function generates the line plots of the phenotype proportions 
    at different distances from a given cell
    '''

    ax.stackplot(dist_bin, cell_propor.T, labels = phenoSet, colors = colors)

    ax.set_xticks(dist_bin)
    # ax.set_xticklabels(['0-25', '26-50', '51-100', '101-150', '151-200'])
    ax.set_xlim([25, 200])
    ax.set_ylim([0, 1])
    ax.set_title(f'Cell #{cell_label}', fontsize = 16)
    ax.set_xlabel('Spatial Bound (\u03BCm)', fontsize = 14)
    ax.set_ylabel('Cell Proportions', fontsize = 14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if legF:
        plt.legend()

    
def plot_mean_neighborhood_profile(ax, selClus, dist_bin, cell_density, phenoSet, maxDens=0.1, legF=0):
    '''This function generates the line plots of the phenotype density 
    at different distances from a given cell
    '''

    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    means = cell_density['means']
    error = cell_density['error']

    if means.ndim == 3:
        maxDens = means[:,:,selClus].max() + error[:,:,selClus].max()

        axesDict = dict()
        for i , key in enumerate(phenoSet):
            plotax = ax.errorbar(dist_bin, means[:, i, selClus], yerr=error[:, i, selClus], color=phenoSet[key])

            ax.set_xticks(dist_bin)
            # ax.set_xticklabels(['0-25', '26-50', '51-100', '101-150', '151-200'])
            ax.set_xlim([20, 210])
            ax.set_ylim([0, maxDens])
            ax.set_title(f'Cluster {selClus}: Densities', fontsize = 16, color = SlTC)
            ax.set_xlabel('Spatial Bound (\u03BCm)', fontsize = 14, color = SlTC)
            ax.set_ylabel('Cell Density', fontsize = 14, color = SlTC)

            ax.set_frame_on(False)
            ax.spines['left'].set_color(SlTC)
            ax.spines['bottom'].set_color(SlTC)
            ax.tick_params(axis='x', colors=SlTC, which='both')
            ax.tick_params(axis='y', colors=SlTC, which='both')
            
            axesDict[key] = plotax

        if legF:
            ax.legend(axesDict.values(), axesDict.keys(),
                    bbox_to_anchor=(-0.05, -0.1), 
                    loc='upper left', 
                    borderaxespad=0, 
                    ncols = 4,
                    facecolor = Sl2BgC,
                    edgecolor = Sl2BgC,
                    labelcolor = SlTC)

def plot_incidence_line(ax, df, phenotype):
    
    SlBgC  = '#0E1117'  # Streamlit Background Color
    SlTC   = '#FAFAFA'  # Streamlit Text Color
    Sl2BgC = '#262730'  # Streamlit Secondary Background Color

    plotax = ax.plot(df.index, df, marker = 'o', markersize = 14, linewidth=2.5)

    ax.set_frame_on(False) # Turn off the Frame
    ax.grid('True', alpha = 0.3)
    ax.spines['left'].set_color(SlTC)
    ax.spines['bottom'].set_color(SlTC)
    ax.tick_params(axis='x', colors=SlTC, which='both')
    ax.tick_params(axis='y', colors=SlTC, which='both')
    ax.legend(plotax, [phenotype], 
              facecolor = Sl2BgC,
              edgecolor = Sl2BgC,
              labelcolor = SlTC,
              fontsize = 16)

def draw_cmp_swatches(color_list):

    swatchFig = plt.figure(figsize = (4,4))
    ax = swatchFig.add_subplot(1,1,1)

    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])

    eL = 1
    buff = 0.2

    row = 0
    for i, color in enumerate(color_list):
        x_pos = (i-row*6)

        ax.add_patch(Rectangle((x_pos+(buff*x_pos), row + buff*row), eL, eL,
                    edgecolor = color,
                    facecolor = color,
                    fill=True))
        
        if ((i+1) % 6) == 0:
            row = row + 1
        
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)