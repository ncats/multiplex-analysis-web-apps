import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import ndimage as ndi


def plot_2d_density(X, Y=None, bins=200, n_pad=40, w=None, ax=None, gaussian_sigma=0.5, cmap=plt.get_cmap('viridis'), vlim=np.array([0.001, 0.98]), circle_type='bg',  box_off=True, return_matrix=False):
    """
    Plot a 2D density plot.

    Parameters:
    - X: array-like, shape (n_samples,)
        The x-coordinates of the data points.
    - Y: array-like, shape (n_samples,), optional
        The y-coordinates of the data points. If None, X is treated as the density matrix.
    - bins: int or array-like, optional
        The number of bins or the bin edges along each dimension.
    - n_pad: int, optional
        The number of padding pixels around the density matrix.
    - w: array-like, shape (n_samples,), optional
        The weights of the data points.
    - ax: matplotlib Axes, optional
        The axes on which to plot the density plot. If None, a new figure and axes will be created.
    - gaussian_sigma: float, optional
        The standard deviation of the Gaussian filter applied to the density matrix.
    - cmap: matplotlib colormap, optional
        The colormap used for the density plot.
    - vlim: array-like, shape (2,), optional
        The lower and upper limits of the colorbar.
    - circle_type: str, optional
        The type of circle to overlay on the density plot. Can be 'bg' (background), 'arch' (arch), or None.
    - box_off: bool, optional
        Whether to turn off the box and ticks of the axes.
    - return_matrix: bool, optional
        Whether to return the density matrix instead of plotting it.

    Returns:
    - If return_matrix is True, returns the density matrix.
    - Otherwise, plots the density matrix on the specified axes.

    """
    if Y is not None:
        if w is not None:
            # Compute the joint histogram and apply Gaussian filter
            b, _, _ = np.histogram2d(X, Y, bins=bins)
            b = ndi.gaussian_filter(b.T, sigma=gaussian_sigma)

            # Compute the weighted histogram and apply Gaussian filter
            s, _, _ = np.histogram2d(X, Y, bins=bins, weights=w)
            s = ndi.gaussian_filter(s.T, sigma=gaussian_sigma)

            # Compute the density matrix
            d = np.zeros_like(b)
            d[b > 0] = s[b > 0] / b[b > 0]
            d = ndi.gaussian_filter(d, sigma=gaussian_sigma)
        else:
            # Compute the histogram and normalize
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
            # Create a circular mask and overlay it on the density matrix
            c = np.meshgrid(np.arange(2 * n_pad + n_bins), np.arange(2 * n_pad + n_bins))
            c = np.sqrt(((c[0] - ((2 * n_pad + n_bins) / 2)) ** 2) + ((c[1] - ((2 * n_pad + n_bins) / 2)) ** 2)) < (0.95 * ((2 * n_pad + n_bins) / 2))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]) + c, vmin=1, vmax=1 + vlim[1], cmap=cmap, shading='gouraud', alpha=1)
        elif circle_type == 'arch':
            # Create an arch-shaped mask and overlay it on the density matrix
            c = (n_bins / 2)
            ax.add_artist(plt.Circle((c + n_pad, c + n_pad), 0.95 * (c + n_pad), color='black', fill=False))
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=-vlim[1], vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)
        else:
            # Plot the density matrix without any overlay
            ax.pcolormesh(np.pad(d, [n_pad, n_pad]), vmin=0, vmax=vlim[1], cmap=cmap, shading='gouraud', alpha=1)

        if box_off is True:
            # Turn off the box and ticks of the axes
            [ax.spines[sp].set_visible(False) for sp in ax.spines]
            ax.set(xticks=[], yticks=[])


def plt_cmap(ax, cmap, extend, width, ylabel):
    """
    Plot a colorbar with a specified colormap.

    Parameters:
    - ax: matplotlib Axes
        The axes on which to plot the colorbar.
    - cmap: matplotlib colormap
        The colormap to use for the colorbar.
    - extend: str
        The extend of the colorbar. Can be 'neither', 'both', 'min', or 'max'.
    - width: float
        The width of the colorbar.
    - ylabel: str
        The label for the colorbar.

    """
    cb = mpl.colorbar.ColorbarBase(ax=ax, cmap=cmap, extend=extend)
    cb.set_ticks([])
    pos = ax.get_position().bounds
    ax.set_position([pos[0], pos[1], width, pos[3]])
    ax.set(ylabel=ylabel)
