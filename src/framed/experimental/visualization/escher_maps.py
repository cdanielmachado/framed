def build_escher_map(fluxes, map_name, abstol=1e-6):
    """ Build escher map for a given flux distribution

    Args:
        fluxes (dict): flux distribution
        map_name (str): name of **escher** map (for a list of maps see *list_escher_maps*)
        abstol (float): tolerance to remove residual fluxes for cleaner visualization (default: 1e-6)

    Returns:
        escher.plots.Builder: escher map object
    """

    import escher

    data = {r_id[2:]: abs(val) if abs(val) > abstol else 0.0 for r_id, val in fluxes.items()}

    colors = [{'type': 'min', 'color': '#f0f0f0', 'size': 8},
              {'type': 'mean', 'color': '#abb7ff', 'size': 20},
              {'type': 'max', 'color': '#0f16fb', 'size': 40}]

    emap = escher.Builder(map_name=map_name, reaction_data=data,
                          reaction_scale=colors,
                          reaction_no_data_size=8,
                          reaction_no_data_color='#ffddda')

    return emap


def list_escher_maps():
    """ List of maps available in **escher**

    Returns:
        list: map names
    """
    import escher

    maps = escher.list_available_maps()
    return [entry['map_name'] for entry in maps]

