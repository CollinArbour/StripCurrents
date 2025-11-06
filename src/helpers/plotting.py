import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker

def add_cms_label(fig, ax=None,
                  cms_text='CMS', status_text='Preliminary',context=None,
                  cms_fontsize=14, status_fontsize=12,context_fontsize=14,
                  y=-1, x_offset=0.01, gap=0.09,
                  top_adjust=0.92):
    """
    Adds a CMS-style label aligned with the plot area (not figure edge).

    Parameters:
        fig             : Matplotlib figure object
        ax              : Optional Axes object; if None, uses plt.gca()
        cms_text        : Bold text (e.g., 'CMS MUON')
        status_text     : Italic text (e.g., 'Preliminary')
        cms_fontsize    : Font size for CMS text
        status_fontsize : Font size for status
        y               : Vertical position in figure coordinates
        x_offset        : Offset *into* the axes region (as fraction of axes width)
        gap             : Horizontal space between CMS and status text
        top_adjust      : Adjusts figure top margin to make room
    """
    if ax is None:
        ax = plt.gca()

    if y < 0:
        y = top_adjust+0.01  # Default to just below the top of the figure

    # Get bounding box of the Axes in figure coordinates
    ax_pos = ax.get_position()
    x_fig = ax_pos.x0 + x_offset * (ax_pos.x1 - ax_pos.x0)

    # Place the CMS and status labels
    fig.text(x_fig, y, cms_text,
             fontsize=cms_fontsize, fontweight='bold',
             ha='left', va='bottom')

    fig.text(x_fig + gap, y, status_text,
             fontsize=status_fontsize, style='italic',
             ha='left', va='bottom')
    
    if context:
        x2_fig = ax_pos.x1 - 0.001 * (ax_pos.x1 - ax_pos.x0)
        # x2_fig = ax_pos.x1 - x_offset * (ax_pos.x1 - ax_pos.x0)
        # fig.text(ax_pos.x1, y, context, ha='right', va='bottom', fontsize=context_fontsize)
        fig.text(x2_fig, y, context, ha='right', va='bottom', fontsize=context_fontsize)

    fig.subplots_adjust(top=top_adjust)

def sci_notation_formatter(x, _):
    if x == 0:
        return "$0$"
    exponent = int(np.floor(np.log10(abs(x))))
    coeff = x / 10**exponent
    return fr"${coeff:.1f} \times 10^{{{exponent}}}$"