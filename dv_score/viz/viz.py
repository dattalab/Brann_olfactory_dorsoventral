import cmcrameri
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats


def dv_cmap(dark=False):
    dv_pal = sns.husl_palette(as_cmap=True)
    if dark:
        hsl = sns.husl_palette(n_colors=256)
        dv_pal = mpl.colors.ListedColormap(hsl[::-1])
    return dv_pal


def dv_colors():
    pal = sns.husl_palette(7)
    return [pal[1], pal[-1]]


def dv_cell_cmap():
    return cmcrameri.cm.romaO


def add_cbar(im, fig, pos=[0.9, 0.2, 0.04, 0.6]):
    cax = fig.add_axes(pos)
    cbar = fig.colorbar(im, cax)
    cbar.solids.set_alpha(1)
    update_cbar(cbar)
    return cbar, cax


def add_horz_cbar(im, fig, pos=[0.2, 0.85, 0.6, 0.05]):
    cb_ax = fig.add_axes(pos)
    cbar = fig.colorbar(im, cax=cb_ax, orientation="horizontal")
    cb_ax.xaxis.set_label_position("top")
    cb_ax.xaxis.set_ticks_position("top")
    cb_ax.xaxis.set_tick_params(pad=0)
    cbar.solids.set_alpha(1)
    update_cbar(cbar)
    return cbar, cb_ax


def update_cbar(cbar, lw=None):
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")
    if lw is None:
        cbar.outline.set_visible(False)
    else:
        cbar.outline.set_linewidth(lw)


def remove_leg_border(legend):
    for lh in legend.legend_handles:
        # none needs to be in quotes
        lh.set_edgecolor("None")


def update_boxen(ax, cls="k", lw=0.5, ls="--", remove_fliers=True):
    """
    Parameters
    ----------
    ax : matplotlib.axes
        Axes containing a boxenplot
    """
    for a in ax.lines:
        a.set_color(cls)
        a.set_linewidth(lw)
        a.set_linestyle(ls)
        a.set_alpha(1)
    for a in ax.collections:
        if isinstance(a, mpl.collections.PatchCollection):
            # remove line surrounding each box
            a.set_linewidth(0)
        else:
            # remove outlier points
            if remove_fliers:
                a.set_alpha(0)


def adjust_lims(
    ax, xmin=None, xmax=None, cls="k", zo=-3, alpha=0.5, add_line=True, **kwargs
):
    vlim = ax.axes.viewLim.extents
    if xmin is None:
        xmin = np.min(vlim[:2])
    if xmax is None:
        xmax = np.max(vlim[2:])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)
    if add_line:
        ax.plot(
            [xmin, xmax],
            [xmin, xmax],
            color=cls,
            ls="--",
            zorder=zo,
            alpha=alpha,
            **kwargs
        )


def offset_pointplot(ax, offset=0.4):
    offset = mpl.transforms.ScaledTranslation(offset, 0, ax.figure.dpi_scale_trans)

    for coll in ax.collections:
        trans = coll.get_transform()
        coll.set_transform(trans + offset)
    # shift everything else:
    for line in ax.lines:
        trans = line.get_transform()
        line.set_transform(trans + offset)


def modify_max(x, thresh=1):
    if x.max() < thresh:
        return 1
    else:
        return x.max()


# for use with plt.step
def get_step(vect, add_max=True):
    vect = np.sort(vect)
    if add_max:
        return (
            np.hstack([vect, modify_max(vect)]),
            np.arange(vect.size + 1) / (vect.size + 1),
        )
    else:
        return vect, np.arange(vect.size) / vect.size


def plot_r(ax, x, y, ha="right", va="bottom", r_loc=(1, 0.02)):
    r, pval = stats.spearmanr(x, y)
    ax.text(
        *r_loc,
        rf"$\rho$:{r:.3f}",
        horizontalalignment=ha,
        verticalalignment=va,
        transform=ax.transAxes,
    )
