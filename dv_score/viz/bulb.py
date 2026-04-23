import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from dv_score.viz import viz


def bulb_colors():
    medial_color = plt.cm.Dark2(0)
    lateral_color = plt.cm.Dark2(3)
    return medial_color, lateral_color


def domains():
    DOMAINS = ["ventro-medial", "dorso-lateral"]
    DOMAIN_NAMES = [d.split("-")[1] for d in DOMAINS]
    new_old_map = dict(zip(DOMAINS, DOMAIN_NAMES))
    return DOMAINS, DOMAIN_NAMES, new_old_map


def plot_c_col(df_geom_dv, plot_kwargs=None, c_col="DV"):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(4.4, 2.0),
        # gridspec_kw={"width_ratios": [ratio, 1]},
        sharey="row",
        sharex=True,
    )
    DOMAINS, _, new_old_map = domains()

    for ax, x in zip(axes, DOMAINS):
        _df = df_geom_dv[df_geom_dv.domain == x]
        im = ax.scatter(_df["z"] / 1e3, _df["y"] / 1e3, c=_df[c_col], **plot_kwargs)
        ax.set_title(f"{new_old_map[x]} domain", y=0.95)
    #     ax.set_xlabel(f"R-C [{MICRON}]")

    # axes[0].set_ylabel(f"D-V [{MICRON}]")
    axes[0].invert_xaxis()
    axes[0].set_xlim(2.65, 0)
    # axes[0].set_aspect("equal")
    sns.despine()
    cax = fig.add_axes([0.5, 0.1, 0.27, 0.04])
    cbar = plt.colorbar(im, cax, orientation="horizontal")
    cbar.set_label(f"{c_col} score (scRNA-seq)", labelpad=4)
    cbar.solids.set_alpha(1)
    cax.xaxis.set_label_position("top")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_tick_params(pad=0)
    cax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(150))
    viz.update_cbar(cbar)

    plt.subplots_adjust(hspace=0.2, wspace=0.15)
    # axes[0].plot([2.5, 2.5], [-1, -2], lw=0.5, color="k")
    # axes[0].plot([2.5, 1.5], [-2, -2], lw=0.5, color="k")
    axes[0].text(2.0, -2.05, "1 mm", ha="center", va="bottom")
    axes[0].text(2.0, -2.2, "A-P axis", ha="center", va="top")
    axes[0].text(2.7, -1.6, "D-V axis", ha="center", va="center", rotation=90)

    start_x, start_y = 2.5, -2.1
    axes[0].arrow(
        start_x,
        start_y,
        0,
        1,
        head_width=0.05,
        head_length=0.05,
        fc="k",
        ec="k",
        linewidth=0.5,
    )
    axes[0].arrow(
        start_x,
        start_y,
        -1,
        0,
        head_width=0.05,
        head_length=0.05,
        fc="k",
        ec="k",
        linewidth=0.5,
    )

    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")

    plt.subplots_adjust(hspace=0, wspace=-0.4)

    return fig, axes
