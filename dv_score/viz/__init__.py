
from importlib import resources as _resources

import matplotlib.pyplot as plt


def moe_style_context():
    style_path = _resources.files(__package__) / "moe_figures.mplstyle"
    return plt.style.context(style_path)