from __future__ import division
import requests
from math import ceil
import numpy as np

from framed.experimental.metanetx import MetaNetX


def iPATH_display(data, output='svg', filename=None, default_opacity=0.1, default_width=2, default_radius=5, **kwargs):

    url = 'https://pathways.embl.de/mapping.cgi'

    elements = edge_mapper(data)

    parameters = {
        'map': 'metabolic',
        'export_type': 'svg',
        'selection': '\n'.join(elements),
        'default_opacity': default_opacity,
        'default_width': default_width,
        'default_radius': default_radius
    }

    parameters.update(kwargs)

    r = requests.post(url, data=parameters)

    if r.headers['Content-Type'] != 'image/svg+xml':
        print(r.text)
        return

    if output == 'svg':
        if filename is None:
            filename = 'iPath.svg'
        with open(filename, 'wb') as f:
            f.write(r.content)

    elif output == 'png':
        from cairosvg import svg2png
        if filename is None:
            filename = 'iPath.png'
        svg2png(bytestring=r.content, write_to=filename)

    elif output == 'notebook':
        from IPython.display import SVG, display_svg
        display_svg(SVG(data=r.content))

    else:
        print('Output format not supported.')


def edge_mapper(data, min_width=0, max_width=25, min_alpha=0.25, max_alpha=1.0):

    if len(data) == 0:
        return 'empty'

    ub = max(data.values())
    lb = min(data.values())

    if ub == lb:
        width_factor = (max_width - min_width) / 2
        alpha_factor = (max_alpha - min_alpha) / 2
        lb = 0
    else:
        width_factor = (max_width - min_width) / (ub - lb)
        alpha_factor = (max_alpha - min_alpha) / (ub - lb)

    edges = []

    for id, value in data.items():
        width = int(ceil((value - lb) * width_factor + min_width))
        alpha = (value - lb) * alpha_factor + min_alpha
        edge = '{:s} W{:d} {:.3f}'.format(id, width, alpha)
        edges.append(edge)

    return edges


def display_fluxes(fluxes, id_convert=None, mnx_path=None, log=False, abstol=1e-6, **kwargs):

    if id_convert is not None:

        if mnx_path is None:
            print('Please specify local path of MetaNetX database.')
            return

        mnx = MetaNetX(mnx_path, version=3)

        # TODO: still need to find how to deal with ambiguous conversions

        fluxes = {
            kegg_id: abs(value)
            for r_id, value in fluxes.items()
            for kegg_id in mnx.translate_reaction_id(r_id[2:], id_convert, 'kegg')
        }
    else:
        fluxes = {key: abs(val) for key, val in fluxes.items()}

    if log:
        lb = np.log10(min([x for x in fluxes.values() if x > abstol]))
        fluxes = {key: np.log10(val) - lb for key, val in fluxes.items() if val > abstol}

    return iPATH_display(fluxes, **kwargs)
