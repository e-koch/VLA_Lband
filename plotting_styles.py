
import matplotlib.pyplot as p
import seaborn as sb

'''
Set style for figures.
'''


def default_figure():
    sb.set(font='Times New Roman', style='ticks')
    sb.set_context("poster")


def twocolumn_figure(fig_ratio=None, font_scale=1.25):
    default_figure()

    width = 8.4
    # Keep the default ratio used in seaborn. This can get overwritten.
    height = (4.4 / 6.4) * width

    figsize = (width, height)

    if fig_ratio is not None:
        figsize = (width, width * fig_ratio)

    sb.set_context("paper", font_scale=font_scale,
                   rc={"figure.figsize": figsize})


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

    sb.set_context("paper", font_scale=font_scale,
                   rc={"figure.figsize": figsize})


def image_shape_ratio(shape):
    '''
    Return the width / height for adjusting figure sizes.
    '''

    return shape[0] / float(shape[1])
