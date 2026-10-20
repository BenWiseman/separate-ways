"""House style for the figures of Separate Ways v6 and its inflation companion.

Computer Modern text and maths, from the TrueType fonts matplotlib ships, so the figures match the
paper's type and embed cleanly in a PDF. Inward ticks on all four sides. One meaning per colour across
every figure: light blue for corrected values and ours, deep purple for published values and the
reference states, orange as a sparing accent for what to look for, greys for guides; text in blue
takes the deep shade. The pairs were chosen by measured colour distance (CIEDE2000) under simulated
deuteranopia, protanopia and tritanopia (Machado et al. 2009) and by lightness, which is what
survives greyscale printing: blue against purple is at worst 40 apart with a lightness gap of 45,
blue against orange at worst 50, and purple against orange at worst 49.
Line styles differ as well as colours.
"""
import matplotlib

BLUE_D, BLUE, BLUE_M, BLUE_L, BLUE_P = "#0b5cad", "#5aaeff", "#2a7fd4", "#c2e0ff", "#eaf4ff"
PURPLE_D, PURPLE, PURPLE_M, PURPLE_L, PURPLE_P = "#2e1159", "#4a1e8c", "#7b5ab5", "#bba9dd", "#eee8f7"
ORANGE_D, ORANGE, ORANGE_L, ORANGE_P = "#b85a0a", "#f28e2b", "#f8bb7c", "#fdebd8"
INK, GREY, GREY_L, SHADE = "#1e1b24", "#6b6874", "#b9b6c0", "#eeecf1"
TEXT_WIDTH = 6.5                                    # inches: the text width of both builds

def apply():
    matplotlib.rcParams.update({
        "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True, "axes.unicode_minus": True,
        "font.size": 9.5, "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "legend.fontsize": 8.2, "legend.frameon": False, "legend.handlelength": 2.4,
        "axes.linewidth": 0.7, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 4.0, "ytick.major.size": 4.0, "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "lines.linewidth": 1.4, "hatch.linewidth": 0.9, "pdf.fonttype": 42, "savefig.dpi": 300})

def panel(ax, label):
    """(a), (b) just above the top-left corner, outside the frame, where no data can sit"""
    ax.text(0.0, 1.025, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5)
