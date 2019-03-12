
import seaborn as sb

'''
Set style for figures.
'''


def default_figure(font_scale=1.0):
    sb.set("poster", font_scale=font_scale,
           palette='colorblind',
           font='Times New Roman', style='ticks',
           rc={'text.usetex': True})


def fullpage_figure(font_scale=1.25):
    default_figure()

    width = 8.4
    height = 11

    figsize = (width, height)

    sb.set("paper", font_scale=font_scale,
           rc={"figure.figsize": figsize,
               'text.usetex': True},
           palette='colorblind',
           font='Times New Roman', style='ticks')


def twocolumn_figure(fig_ratio=None, font_scale=1.25):
    default_figure()

    width = 8.4
    # Keep the default ratio used in seaborn. This can get overwritten.
    height = (4.4 / 6.4) * width

    figsize = (width, height)

    if fig_ratio is not None:
        figsize = (width, width * fig_ratio)

    sb.set("paper", font_scale=font_scale,
           rc={"figure.figsize": figsize,
               'text.usetex': True},
           palette='colorblind',
           font='Times New Roman', style='ticks')


def twocolumn_twopanel_figure(fig_ratio=None, **kwargs):

    # Use half the ratio of a one column figure.
    if fig_ratio is None:
        fig_ratio = (4.4 / 6.4) / 2

    twocolumn_figure(fig_ratio=fig_ratio, **kwargs)


def onecolumn_figure(fig_ratio=None, font_scale=1.2):
    '''
    fig_ratio is width / height.
    '''
    default_figure()

    # About the width (in inches) of a column
    width = 4.2
    # Keep the default ratio used in seaborn. This can get overwritten.
    height = (4.4 / 6.4) * width

    figsize = (width, height)

    if fig_ratio is not None:
        figsize = (width, width * fig_ratio)

    sb.set("paper", font_scale=font_scale,
           rc={"figure.figsize": figsize,
               'text.usetex': True},
           palette='colorblind',
           font='Times New Roman', style='ticks')
    # sb.set_palette("colorblind")


def onecolumn_Npanel_figure(N, font_scale=1.2):

    width = 4.2

    height = (4.4 / 6.4) * N * width

    onecolumn_figure(fig_ratio=height / width, font_scale=font_scale)


def onecolumn_twopanel_figure(font_scale=1.2):

    onecolumn_Npanel_figure(N=2, font_scale=font_scale)


def image_shape_ratio(shape):
    '''
    Return the width / height for adjusting figure sizes.
    '''

    return shape[0] / float(shape[1])


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2, (y1 - y2) / 2, v2)
    adjust_yaxis(ax1, (y2 - y1) / 2, v1)


def adjust_yaxis(ax, ydif, v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny > maxy or (-miny == maxy and dy > 0):
        nminy = miny
        nmaxy = miny * (maxy + dy) / (miny + dy)
    else:
        nmaxy = maxy
        nminy = maxy * (miny + dy) / (maxy + dy)
    ax.set_ylim(nminy + v, nmaxy + v)
