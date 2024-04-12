# Author: Alex Baras, MD, PhD (https://github.com/alexbaras)
# NCATS Maintainer: Dante J Smith, PhD (https://github.com/djsmith17)
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
import seaborn as sns
import altair as alt
alt.data_transformers.disable_max_rows()
from scipy import ndimage as ndi
from scipy.stats import binned_statistic_2d

def plot_2d_density(X, Y=None, bins=200, n_pad=40, w=None, ax=None, gaussian_sigma=0.5, cmap=plt.get_cmap('viridis'), vlim=np.array([0.001, 0.98]), circle_type='bg', box_off=True, return_matrix=False, legendtype = 'colorbar'):
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

    if np.isscalar(bins):
        n_bins = bins
    else:
        n_bins = len(bins[0]) - 1

    if Y is not None:
        if w is not None:
            b, _, _ = np.histogram2d(X, Y, bins=bins)
            b = ndi.gaussian_filter(b.T, sigma=gaussian_sigma)

            s, xedges, yedges = np.histogram2d(X, Y, bins=bins, weights=w)
            s = ndi.gaussian_filter(s.T, sigma=gaussian_sigma)

            d = np.zeros_like(b)
            # d[b > 0] = s[b > 0] / b[b > 0]
            d = s
            d = ndi.gaussian_filter(d, sigma=gaussian_sigma)
        else:
            d, xedges, yedges = np.histogram2d(X, Y, bins=bins)
            d /= np.sum(d)
            d = ndi.gaussian_filter(d.T, sigma=gaussian_sigma)
    else:
        d = X

    if return_matrix:
        x_bin_indices = np.digitize(X, xedges[:-1])-1
        y_bin_indices = np.digitize(Y, yedges[:-1])-1
        # bin_indices = [(indx*n_bins + indy) for (indx, indy) in zip(x_bin_indices, y_bin_indices)]

        # bin_indices_df = pd.DataFrame(bin_indices, columns = ['bin_num'])
        # bin_indices_df['index'] = bin_indices_df.index
        # bin_indices_df_group = bin_indices_df.groupby('bin_num')['index'].apply(list)

        # tuple_list = [(indx, indy) for (indx, indy) in zip(x_bin_indices, y_bin_indices)]

        bin_indices_df = pd.DataFrame(data = {'indx': x_bin_indices.flatten(),
                                              'indy': y_bin_indices.flatten(),
                                              'valx': X,
                                              'valy': Y})
        return d, bin_indices_df
    else:
        if d[d > 0].shape == (0,):
            vmin = 0
            vmax = 1
            vlim = [vmin, vmax]
        else:
            if ~all((d < 0)[0]):
                if np.isscalar(vlim):
                    vlim = np.array([0, np.quantile(d[d > 0].flatten(), vlim)])
                else:
                    if np.all((vlim < 1) & (vlim > 0)):
                        vlim = np.quantile(d[d > 0].flatten(), vlim)
            else:
                vlim = np.zeros(200)

        if ax is None:
            _, ax = plt.subplots()

        if np.isscalar(bins):
            n_bins = bins
        else:
            n_bins = len(bins[0]) - 1

        if circle_type == 'bg':
            extend = 'max'
            c = np.meshgrid(np.arange(2 * n_pad + n_bins), np.arange(2 * n_pad + n_bins))
            c = np.sqrt(((c[0] - ((2 * n_pad + n_bins) / 2)) ** 2) + ((c[1] - ((2 * n_pad + n_bins) / 2)) ** 2)) < (0.95 * ((2 * n_pad + n_bins) / 2))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]) + c, vmin=1, vmax=1 + vlim[1], cmap=cmap, shading='gouraud', alpha=1)
            # ax.pcolormesh(np.log10(np.pad(d, [n_pad, n_pad]) + c + 1), vmin=np.log10(2), vmax=np.log10(2 + vlim[1]), cmap=cmap, shading='gouraud', alpha=1)
        elif circle_type == 'arch':
            extend = 'both'
            # if any((d<0)[0]):
            #     vmin = np.quantile(d[d < 0].flatten(), 0.03)
            # else:
            #     vmin = 0.03
            # if any((d>0)[0]):
            #     vmax = np.quantile(d[d > 0].flatten(), 0.97)
            # else:
            #     vmax = 0.97
            c = (n_bins / 2)
            ax.add_artist(plt.Circle((c + n_pad, c + n_pad), 0.95 * (c + n_pad), color='black', fill=False))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=-vlim[1], vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)
        else:
            extend = 'max'
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=0, vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)

        # Create the color bar
        if legendtype == 'colorbar':
            cax = ax.inset_axes([0.95, 0.1, 0.01, 0.85])
            cmap_lim = None # [np.min(d), np.max(d)]
            plt_cmap(ax=cax, cmap=cmap, extend=extend, width=0.01, lim = cmap_lim)
        elif legendtype == 'legend':
            cax = ax.inset_axes([0.95, 0.1, 0.01, 0.85])
            cmap_lim = [-3, -2, -1, 0, 1, 2, 3]
            plt_cmap(ax=cax, cmap=cmap, extend=extend, width=0.01, lim = cmap_lim)


        if box_off is True:
            [ax.spines[sp].set_visible(False) for sp in ax.spines]
            ax.set(xticks=[], yticks=[])


def plt_cmap(ax, cmap, extend, width, lim = None, ylabel = None):
    '''
    plt_cmap(ax, cmap, extend, width, ylabel) draws a colorbar 
    for the current colormap at the correct
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
    cb = mpl.colorbar.Colorbar(ax=ax, cmap=cmap, extend=extend)
    cb.set_ticks([])
    pos = ax.get_position().bounds
    ax.set_position([pos[0], pos[1], width, pos[3]])

    if lim is not None:
        cb.set_ticks([0, 1])
        cb.set_ticklabels([lim[0], lim[1]])

    if ylabel is not None:
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
    '''
    This function generates the line plots of the phenotype density 
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

def plot_mean_neighborhood_profile(ax, dist_bin, npf_dens_mean, cluster_title, max_dens=0.1, leg_flag=0):
    '''
    This function generates the line plots of the phenotype density 
    at different distances from a given cell
    '''

    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    tab20 = plt.get_cmap('tab20')
    tab20_new = ListedColormap(tab20(np.arange(256)))

    phenotypes = npf_dens_mean['phenotype'].unique()
    axesDict = dict()
    for ii, phenotype in enumerate(phenotypes):
        # Find the phenotype in the dataframe
        npf_dens_mean_pheno = npf_dens_mean[npf_dens_mean['phenotype'] == phenotype]

        plotax = ax.errorbar(x = npf_dens_mean_pheno.dist_bin,
                             y = npf_dens_mean_pheno.density_mean,
                             yerr=npf_dens_mean_pheno.density_sem,
                             color=tab20_new(ii))
        axesDict[phenotype] = plotax

    plt.axhline(y=0, color='w', linestyle='--')

    ax.set_xticks(dist_bin)
    ax.set_xlim([0, 225])
    ax.set_ylim(max_dens)
    ax.set_title(cluster_title, fontsize = 20, color = slc_text)
    ax.set_xlabel('Spatial Bound (\u03BCm)', fontsize = 14, color = slc_text)
    ax.set_ylabel('Cell Density', fontsize = 14, color = slc_text)

    # ax.set_frame_on(False)
    ax.spines[['left', 'bottom']].set_color(slc_text)
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='x', colors=slc_text, which='both')
    ax.tick_params(axis='y', colors=slc_text, which='both')

    if leg_flag:
        ax.legend(axesDict.values(), axesDict.keys(),
                  bbox_to_anchor=(-0.05, -0.1),
                  loc='upper left',
                  fontsize = 12,
                  borderaxespad=0,
                  ncols = 4,
                  facecolor = slc_bg,
                  edgecolor = slc_bg,
                  labelcolor = slc_text)

def plot_incidence_line(ax, df, phenotype):
    '''
    Plot the incidence of a given phenotype over time

    Args:
        ax: Matplotlib axis handle
        df: Pandas DataFrame
        phenotype: String of the phenotype to plot

    Returns:
        None    
    '''

    # Streamlit Theming
    slc_bg   = '#0E1117'  # Streamlit Background Color
    slc_text = '#FAFAFA'  # Streamlit Text Color
    slc_bg2  = '#262730'  # Streamlit Secondary Background Color

    plotax = ax.plot(df.index, df, marker = 'o', markersize = 14, linewidth=2.5)

    ax.set_frame_on(False) # Turn off the Frame
    ax.grid('True', alpha = 0.3)
    ax.spines['left'].set_color(slc_text)
    ax.spines['bottom'].set_color(slc_text)
    ax.tick_params(axis='x', colors=slc_text, which='both')
    ax.tick_params(axis='y', colors=slc_text, which='both')
    ax.legend(plotax, [phenotype],
              facecolor = slc_bg2,
              edgecolor = slc_bg2,
              labelcolor = slc_text,
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