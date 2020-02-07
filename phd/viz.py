import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.path import Path
from matplotlib.patches import BoxStyle
from matplotlib.offsetbox import AnchoredText
import matplotlib
import altair as alt
import bokeh.themes
import bokeh.io
import altair as alt


def color_palette():
    """ 
    Returns my preferred color scheme 
    """
    colors = {
        "dark_purple": "#5F2E88",
        "purple": "#7E59A2",
        "light_purple": "#A17DB8",
        "pale_purple": "#dfd6e5",
        "dark_orange": "#F38227",
        "orange": "#E39943",
        "light_orange": "#EEBA7F",
        "pale_orange": "#f2d4b6",
        "dark_blue": "#3F60AC",
        "blue": "#7292C7",
        "light_blue": "#A5B3CC",
        "pale_blue": "#dae4f1",
        "dark_red": "#9C372F",
        "red": "#C76A6A",
        "light_red": "#E39C9D",
        "pale_red": "#edcccc",
        "dark_green": "#395A34",
        "green": "#688A2F",
        "light_green": "#B3CD86",
        "pale_green": "#d8e2c3",
        "dark_brown": "#764f2a",
        "brown": "#c2996c",
        "light_brown": "#e1bb96",
        "pale_brown": "#efccaf",
        "black": "#444147",
        "grey": "#EFEFEF",
        "gray": "#EFEFEF",
        "light_grey": "#6D6F72",
        "light_gray": "#6D6F72",
    }
    palette = [
        v
        for k, v in colors.items()
        if (v not in ["grey", "gray", "dark_purple", "light_grey"])
        and ("pale" not in [v])
    ]
    return colors, palette


def phd_style():
    """
    Sets the plotting style to my preference
    """
    rc = {
        "axes.facecolor": "#EFEFEF",
        "font.family": "sans-serif",
        "font.family": "Myriad Pro",
        "font.style": "normal",
        "pdf.fonttype": 42,
        "axes.edgecolor": "#444147",
        "axes.labelcolor": "#444147",
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.axisbelow": True,
        "axes.linewidth": 0.25,
        "axes.titlesize": 8,
        "text.color": "#444147",
        "axes.grid": True,
        "lines.linewidth": 0.75,
        "lines.dash_capstyle": "round",
        "grid.linestyle": "-",
        "grid.linewidth": 0.75,
        "grid.color": "#ffffff",
        "axes.labelsize": 8,
        "xtick.color": "#444147",
        "ytick.color": "#444147",
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "xtick.major.pad": 6,
        "ytick.major.pad": 6,
        "xtick.minor.size": 0,
        "ytick.minor.size": 0,
        "legend.fontsize": 6,
        "legend.frameon": True,
        "axes.xmargin": 0.01,
        "axes.ymargin": 0.01,
        "figure.facecolor": "white",
        "figure.dpi": 300,
        "errorbar.capsize": 1,
        "savefig.bbox": "tight",
    }

    plt.rc("text.latex", preamble=r"\usepackage{mathpazo}")
    matplotlib.style.use(rc)
    colors, palette = color_palette()
    sns.set_palette(palette)
    return colors, palette


def bokeh_theme(return_color_list=True):
    theme_json = {
        "attrs": {
            "Figure": {"background_fill_color": "#EEEEEE",},
            "Axis": {
                "axis_line_color": "slategray",
                "major_tick_line_color": None,
                "minor_tick_line_color": None,
            },
            "Legend": {
                "border_line_color": "slategray",
                "background_fill_color": "#EEEEEE",
                "border_line_width": 0.75,
                "background_fill_alpha": 0.75,
            },
            "Grid": {"grid_line_color": "#FFFFFF", "grid_line_width": 0.75,},
            "Text": {
                "text_font_style": "regular",
                "text_font": "Myriad Pro",
                "text_font_size": 10,
            },
            "Title": {
                "background_fill_color": "#EEEEEE",
                "text_color": "#3c3c3c",
                "align": "center",
                "text_font": "Myriad Pro",
                "offset": 3,
            },
        }
    }

    colors, palette = color_palette()
    theme = bokeh.themes.Theme(json=theme_json)
    bokeh.io.curdoc().theme = theme
    if return_color_list:
        return [colors, palette]
    else:
        return colors

def altair_theme():
    """
    Sets a theme for the plotting library Altair to match the style of my PhD.
    """

    colors, palette = color_palette()
    def _theme():
        return {
            'config': {
                'background': 'white',
                    'group': { 
                    'fill': 'white', 
                    },
                'view': {
                    'strokeWidth': 0,
                    'height': 300,
                    'width': 400,
                    'fill': '#EEEEEE'
                    },
                'mark': {
                    'strokeWidth': 0.5,
                    'stroke': 'black'
                },
                'axis': {
                    'domainColor': colors['black'],
                    'domainWidth': 0.5,
                    'labelColor': colors['black'],
                    'labelFontSize': 8,
                    'labelFont': 'Myriad Pro',
                    'titleFont': 'Myriad Pro',
                    'titleFontWeight': 400,
                    'titleFontSize': 8,
                    'grid': True,
                    'gridColor': 'white',
                    'gridWidth': 0.75,
                    'ticks': True,
                    'tickColor': 'white',
                    'tickOffset': 8,
                    'tickWidth': 1.5
                },
                'range': {
                    'category': palette
                },
                'legend': {
                    'labelFontSize': 8,
                    'labelFont': 'Myriad Pro',
                    'titleFont': 'Myriad Pro',
                    'titleFontSize': 8,
                    'titleFontWeight': 400
                },
                'title' : { 
                    'font': 'Myriad Pro',
                    'fontWeight': 400,
                    'fontSize': 8,
                    'anchor': 'middle'
                }
                  }
                }

    alt.themes.register('phd', _theme)# enable the newly registered theme
    alt.themes.enable('phd')
    return colors, palette


def color_selector(style):
    """
    Select the color palette of your choice.

    Parameters
    ----------
    style: str "mut" or "pboc"
        A string identifier for the style. "mut" gives colors for single and double mutants.
        "pboc" returns the PBoC2e color palette.

    Returns
    -------
    colors: dict
        Dictionary of colors. If "dna", "double", or "inducer" is the selected style,
        keys will be the mutants in upper case. Double mutant keys will be DNA-IND. For
        pboc, the keys will be the typical color descriptors. 

    """
    # Ensure the provided style name makes sense.
    if style.lower() not in ["mut", "pboc"]:
        raise ValueError(
            "Provided style must be 'pboc' or 'mut'. {} provided.".format(style)
        )

    # Set the color styles and return.
    if style.lower() == "mut":
        colors = {
            "Y20I": "#738FC1",
            "Q21A": "#7AA974",
            "Q21M": "#AB85AC",
            "F164T": "#A97C50",
            "Q294K": "#5D737E",
            "Q294V": "#D56C55",
            "Q294R": "#B2AF58",
            "Y20I-F164T": "#2d98da",
            "Y20I-Q294K": "#34495e",
            "Y20I-Q294V": "#8854d0",
            "Q21A-F164T": "#4b6584",
            "Q21A-Q294K": "#EE5A24",
            "Q21A-Q294V": "#009432",
            "Q21M-F164T": "#1289A7",
            "Q21M-Q294K": "#6F1E51",
            "Q21M-Q294V": "#006266",
            "WT": "#3C3C3C",
        }

    elif style.lower() == "pboc":
        colors = {
            "green": "#7AA974",
            "light_green": "#BFD598",
            "pale_green": "#DCECCB",
            "yellow": "#EAC264",
            "light_yellow": "#F3DAA9",
            "pale_yellow": "#FFEDCE",
            "blue": "#738FC1",
            "light_blue": "#A9BFE3",
            "pale_blue": "#C9D7EE",
            "red": "#D56C55",
            "light_red": "#E8B19D",
            "pale_red": "#F1D4C9",
            "purple": "#AB85AC",
            "light_purple": "#D4C2D9",
            "dark_green": "#7E9D90",
            "dark_brown": "#905426",
        }
    return colors


def titlebox(
    ax, text, color, bgcolor=None, size=8, boxsize="10%", pad=0.02, loc=10, **kwargs
):
    """Sets a colored box about the title with the width of the plot"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size=boxsize, pad=pad)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.spines["top"].set_visible(True)
    cax.spines["right"].set_visible(True)
    plt.setp(cax.spines.values(), color=color)
    if bgcolor != None:
        cax.set_facecolor(bgcolor)
    else:
        cax.set_facecolor("white")
    at = AnchoredText(text, loc=loc, frameon=False, prop=dict(size=size, color=color))
    cax.add_artist(at)


def ylabelbox(ax, text, color, bgcolor=None, size=6, boxsize="15%", pad=0.02, **kwargs):
    """Sets a colored box about the title with the width of the plot"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("left", size=boxsize, pad=pad)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.spines["top"].set_visible(True)
    cax.spines["right"].set_visible(True)
    plt.setp(cax.spines.values(), color=color)
    if bgcolor != None:
        cax.set_facecolor(bgcolor)
    else:
        cax.set_facecolor("white")

    at = AnchoredText(
        text,
        loc=10,
        frameon=False,
        prop=dict(rotation="vertical", size=size, color=color),
    )
    cax.add_artist(at)
