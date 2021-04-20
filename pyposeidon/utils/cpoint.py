from scipy.spatial.distance import cdist


def closest_node(node, nodes):
    return nodes[cdist([node], nodes).argmin()]
