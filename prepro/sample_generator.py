#
#
#

import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d, HalfspaceIntersection, ConvexHull
from shapely.geometry import Polygon, box, Point, LineString
from scipy.spatial.distance import cdist
import math
import sys


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite Voronoi regions in a 2D diagram to finite regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'
    Returns
    -------
    regions : list of tuples
        Regions as list of tuples of (x,y) pairs
    """

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        pt_min = vor.points.min(axis=0)
        pt_max = vor.points.max(axis=0)
        radius = np.linalg.norm(pt_max - pt_min) * 10.0

    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 not in region:
            new_regions.append([vor.vertices[i].tolist() for i in region])
            continue

        ridges = all_ridges[p_idx]
        pts = []
        for p2, v1, v2 in ridges:
            if v2 < 0 or v1 < 0:
                v_finite = v1 if v1 >= 0 else v2
                finite_vertex = vor.vertices[v_finite]
                tangent = vor.points[p2] - vor.points[p_idx]
                tangent = tangent / np.linalg.norm(tangent)
                normal = np.array([-tangent[1], tangent[0]])

                direction = np.sign(np.dot(finite_vertex - center, normal)) * normal
                far_point = finite_vertex + direction * radius
                pts.append(tuple(finite_vertex))
                pts.append(tuple(far_point))
            else:
                pts.append(tuple(vor.vertices[v1]))
                pts.append(tuple(vor.vertices[v2]))

        pts_unique = np.array(list({p for p in pts}))
        if pts_unique.shape[0] < 3:
            new_regions.append(pts_unique.tolist())
            continue
        angles = np.arctan2(
            pts_unique[:, 1] - vor.points[p_idx, 1],
            pts_unique[:, 0] - vor.points[p_idx, 0],
        )
        order = np.argsort(angles)
        poly_pts = pts_unique[order].tolist()
        new_regions.append(poly_pts)

    return new_regions


def bray_curtis_distance(a, b):
    """
    Compute the Bray-Curtis distance between two arrays of points.

    The Bray-Curtis distance is a measure of dissimilarity between two
    probability distributions. It is defined as the sum of the absolute
    differences between the two distributions, divided by the sum of the
    absolute values of the two distributions.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points

    Returns
    -------
    distances : array
        Array of distances between points in `a` and `b`
    """
    return np.sum(np.abs(a - b), axis=1) / np.sum(np.abs(a + b), axis=1)


def directional_distance(a, b, x_weight=1.0, y_weight=1.0):
    """
    Compute the Euclidean distance between two arrays of points, but with
    different weights for the x and y coordinates.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points
    x_weight, y_weight : float
        Weights for the x and y coordinates, respectively

    Returns
    -------
    distances : array
        Array of distances between points in `a` and `b`
    """
    dx = np.abs(a[:, 0] - b[0]) * x_weight
    dy = np.abs(a[:, 1] - b[1]) * y_weight
    return np.sqrt(dx**2 + dy**2)


def canberra_distance(a, b):
    """
    Compute the Canberra distance between two arrays of points.

    The Canberra distance is a measure of dissimilarity between two
    points in a multi-dimensional space. It is defined as the sum of
    the absolute differences between the two points, divided by the sum
    of the absolute values of the two points.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points

    Returns
    -------
    distances : array
        Array of distances between points in `a` and `b`
    """
    numerator = np.abs(a - b)
    denominator = np.abs(a) + np.abs(b)
    return np.sum(np.where(denominator != 0, numerator / denominator, 0), axis=1)


def mahalanobis_distance(a, b, cov_inv):
    """
    Compute the Mahalanobis distance between two arrays of points.

    The Mahalanobis distance is a measure of dissimilarity between two
    points in a multi-dimensional space. It is defined as the square root
    of the sum of the squared differences between the two points, divided
    by the covariance matrix of the points.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points
    cov_inv : array-like
        Inverse of the covariance matrix of the points

    Returns
    -------
    distances : array
        Array of distances between points in `a` and `b`
    """
    delta = a - b
    return np.sqrt(np.sum(delta @ cov_inv * delta, axis=1))


def bregman_distance(a, b, convex_func, convex_grad):
    """
    Compute the Bregman distance between two arrays of points.

    The Bregman distance is a measure of dissimilarity between two points
    in a multi-dimensional space. It is defined as the difference between
    the value of a convex function at the two points, minus the dot product
    of the difference between the two points and the gradient of the convex
    function at the second point.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points
    convex_func : callable
        Convex function to compute the Bregman distance
    convex_grad : callable
        Gradient of the convex function

    Returns
    -------
    distances : array
        Array of distances between points in `a` and `b`
    """
    F_diff = convex_func(a) - convex_func(b)
    linear_term = np.sum((a - b) * convex_grad(b), axis=1)
    return F_diff - linear_term


def x_distance(a, b):
    """
    Compute the absolute difference between the x-coordinates of two arrays of points.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points

    Returns
    -------
    distances : array
        Array of absolute differences between the x-coordinates of points in `a` and `b`
    """
    return np.abs(a[:, 0] - b[0])


def y_distance(a, b):
    """
    Compute the absolute difference between the y-coordinates of two arrays of points.

    Parameters
    ----------
    a, b : array-like
        Input arrays of points

    Returns
    -------
    distances : array
        Array of absolute differences between the y-coordinates of points in `a` and `b`
    """
    return np.abs(a[:, 1] - b[1])


def bounded_voronoi(points, bounds):
    """
    Compute the bounded Voronoi regions of a set of points within a given box.

    Parameters
    ----------
    points : array-like
        Input array of points
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box

    Returns
    -------
    regions : list
        List of bounded Voronoi regions. If a region is infinite, it is represented
        as None.
    """
    vor = Voronoi(points)
    min_x, min_y, max_x, max_y = bounds
    box_poly = box(min_x, min_y, max_x, max_y)

    regions = []
    for region_index in vor.point_region:
        region = vor.regions[region_index]
        if not region or -1 in region:  # skip infinite
            regions.append(None)
            continue
        polygon = Polygon([vor.vertices[i] for i in region])
        poly_clipped = polygon.intersection(box_poly)
        if not poly_clipped.is_empty:
            regions.append(poly_clipped)
        else:
            regions.append(None)
    return regions


def laguerre_voronoi(points, weights, bounds):
    """
    Compute the Laguerre Voronoi regions of a set of points within a given box.

    Parameters
    ----------
    points : array-like
        Input array of points
    weights : array-like
        Input array of weights associated to each point
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box

    Returns
    -------
    regions : list
        List of bounded Laguerre Voronoi regions. If a region is infinite, it is represented
        as None.
    """
    pts = np.asarray(points)
    w = np.asarray(weights)
    min_x, min_y, max_x, max_y = bounds
    S = np.sum(pts**2, axis=1) - w
    regions = []

    from shapely.geometry import box

    box_poly = box(min_x, min_y, max_x, max_y)

    for i, (xi, yi) in enumerate(pts):
        hs = []

        hs += [[-1, 0, min_x], [1, 0, -max_x], [0, -1, min_y], [0, 1, -max_y]]

        Si = xi * xi + yi * yi - w[i]

        for j, (xj, yj) in enumerate(pts):
            if i == j:
                continue
            a, b = 2 * (xj - xi), 2 * (yj - yi)
            c = S[j] - Si
            hs.append([a, b, -c])

        halfspaces = np.array(hs)

        interior = np.array([xi, yi])
        try:
            hs_int = HalfspaceIntersection(halfspaces, interior)
            verts = hs_int.intersections
            if len(verts) == 0:
                regions.append(None)
                continue
            hull = ConvexHull(verts)
            poly = Polygon(verts[hull.vertices])
            regions.append(poly)
        except Exception:
            regions.append(None)

    return regions


def centroid(poly):
    """
    Compute the centroid of a shapely Polygon.

    Parameters
    ----------
    poly : shapely Polygon
        Input polygon

    Returns
    -------
    centroid : array-like
        (x, y) coordinates of the centroid of the polygon if it is finite and non-empty.
        Otherwise, return None.
    """

    if poly is None or poly.area == 0:
        return None
    c = poly.centroid
    return np.array([c.x, c.y])


def step_lloyd(
    points,
    bounds,
    bundles,
    alpha=1.0,
    distance_metric="euclidean",
    beta_outside=0.0,
    finite_radius=None,
):
    """
    One step of the Lloyd algorithm.

    Parameters
    ----------
    points : array-like
        Array of shape (n, 2) of points to be moved.
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box.
    bundles : array-like
        Array of shape (m, 2) of bundle points.
    alpha : float
        Weight of the bundle points.
    distance_metric : str, optional
        Distance metric to use. Currently only 'euclidean' and 'mahalanobis' are supported.
    beta_outside : float, optional
        Weight of the points outside the domain.
    finite_radius : float, optional
        Radius to used for the finite Voronoi regions.

    Returns
    -------
    new_points : array-like
        Array of shape (n, 2) of the new positions of the points.
    """
    regions_in_bounds = bounded_voronoi(points, bounds)
    new_points = []

    vor = Voronoi(points)
    finite_polys_pts = voronoi_finite_polygons_2d(vor, radius=finite_radius)

    min_x, min_y, max_x, max_y = bounds
    domain_poly = box(min_x, min_y, max_x, max_y)

    cov_matrix = np.array([[1.0, 0.5], [0.0, 100]])
    if distance_metric == "mahalanobis":
        cov_inv = np.linalg.inv(cov_matrix)

    finite_shapely = []
    for pts in finite_polys_pts:
        try:
            poly = Polygon(pts)
            if not poly.is_valid:
                poly = poly.buffer(0)
            finite_shapely.append(poly)
        except Exception:
            finite_shapely.append(None)

    for i, region_in in enumerate(regions_in_bounds):
        # centroïde de la partie à l'intérieur (bounded_voronoi renvoie la partie dans D)
        c_tmp = centroid(region_in)
        if c_tmp is None:
            c_inside = points[i].copy()
        else:
            c_inside = np.asarray(c_tmp)

        if distance_metric == "euclidean":
            dist_to_bundles = np.linalg.norm(bundles - c_inside, axis=1)
        elif distance_metric == "mahalanobis":
            dist_to_bundles = mahalanobis_distance(
                bundles, np.array([c_inside]), cov_inv
            )
        else:
            # Y : J'ai enlevé les autres normes pour le moment ; Frank, tu pourras les rajouter si tu veux
            dist_to_bundles = np.linalg.norm(bundles - c_inside, axis=1)

        closest_bundle = bundles[np.argmin(dist_to_bundles)]

        # Y : Truc important ici : calcul des aires intérieures et extérieures
        poly_total = finite_shapely[i]
        if poly_total is None or poly_total.is_empty:
            # fallback on va dire raisonnable : considérer juste la partie dans D (on ne pénalise pas trop)
            area_total = region_in.area if region_in is not None else 0.0
            area_inside = region_in.area if region_in is not None else 0.0
        else:
            # clipper le polygone "total" avec une grosse boîte pour stabiliser (déjà fait dans voronoi_finite_polygons_2d)
            # aire totale est l'aire de poly_total (dans la grande boîte)
            area_total = float(poly_total.area)
            # aire à l'intérieur du domaine = intersection poly_total ∩ domain_poly
            inter = poly_total.intersection(domain_poly)
            area_inside = (
                float(inter.area) if (inter is not None and not inter.is_empty) else 0.0
            )

        area_out = max(area_total - area_inside, 0.0)

        # normalisation optionnelle pour garder sens des coefficients (utile si tu changes l'échelle du domaine)
        # area_out = area_out / (domain_poly.area)   #Y : décommente pas pour l'instant

        beta_local = beta_outside * area_out

        if beta_local > 0:
            new_pt = (c_inside + alpha * closest_bundle + beta_local * c_inside) / (
                1.0 + alpha + beta_local
            )
        else:
            new_pt = (c_inside + alpha * closest_bundle) / (1.0 + alpha)

        new_points.append(new_pt)

    return np.array(new_points)


def plot_state(points, regions, bundles, bounds):
    """
    Plot the current state of the Lloyd algorithm with the given points, regions, bundles, and bounds.

    Parameters
    ----------
    points : array-like
        The points to plot
    regions : list of shapely Polygon
        The Voronoi regions corresponding to the points, clipped to the bounds
    bundles : array-like
        The locations of the bundles
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box
    """
    fig, ax = plt.subplots(facecolor="black")
    ax.set_facecolor("black")

    for region in regions:
        if region and region.is_valid and not region.is_empty:
            x, y = region.exterior.xy
            ax.plot(x, y, color="white", linewidth=1.0)

    ax.set_xlim(0.25 + bounds[0], bounds[2] - 0.25)
    ax.set_ylim(0.25 + bounds[1], bounds[3] - 0.25)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.show()


def option(A, B, epsilon=0.5):
    """
    Fusionne les points A et B en regroupant les points de A qui sont à une distance inférieure à epsilon d'un point de B en un seul point.

    Parameters
    ----------
    A : array-like
        Les points à fusionner
    B : array-like
        Les points de référence
    epsilon : float, optional
        La distance maximale pour fusionner un point de A avec un point de B

    Returns
    -------
    A_fusionne : array-like
        Les points de A fusionnés avec les points de B
    """

    clustering = DBSCAN(eps=epsilon, min_samples=1).fit(A)
    labels = clustering.labels_
    A_fusionne = []

    for label in np.unique(labels):
        cluster = A[labels == label]

        dists = cdist(cluster, B)
        min_dist = np.min(dists)

        if min_dist < epsilon:
            barycentre = cluster.mean(axis=0)
            A_fusionne.append(barycentre)
        else:
            A_fusionne.extend(cluster)

    return np.array(A_fusionne)


def option_2(A, B, epsilon=0.5):
    """
    Fusionne les points A et B en remplaçant les points de A qui sont à une distance inférieure à epsilon d'un point de B par ce point de B.
    Ensuite, ajoute des points supplémentaires à une distance epsilon alentour de chaque point de B, dans des directions équidistantes.

    Parameters
    ----------
    A : array-like
        Les points à fusionner
    B : array-like
        Les points de référence
    epsilon : float, optional
        La distance maximale pour fusionner un point de A avec un point de B

    Returns
    -------
    A_fusionne : array-like
        Les points de A fusionnés avec les points de B et les nouveaux points alentour
    """

    A_new = A.copy()

    distances = cdist(A, B)  # shape: (len(A), len(B))

    for i, dists in enumerate(distances):
        # Trouver l'indice du bundle le plus proche
        j_min = np.argmin(dists)
        if dists[j_min] < epsilon:
            A_new[i] = B[j_min]  # Remplacer par le point du bundle

    A_unique = np.unique(A_new, axis=0)

    n_directions = 8
    nouveaux_points = []

    for point in B:
        angles = np.linspace(0, 2 * np.pi, n_directions, endpoint=False)
        rayon = epsilon + 1e-4  # epsilon + 0.0001
        for theta in angles:
            x = point[0] + rayon * np.cos(theta)
            y = point[1] + rayon * np.sin(theta)
            nouveaux_points.append([x, y])

    points_autour = np.array(nouveaux_points)

    A_final = np.vstack((A_unique, points_autour))

    return A_final


def run_lloyd(
    num_points=10000,
    num_bundle=50,
    alpha=0.015,
    iterations=20,
    distance_metric="euclidean",
):  # 10000 pour poplar
    """
    Run the Lloyd algorithm to create a point distribution.

    Parameters
    ----------
    num_points : int
        Number of points in the distribution.
    num_bundle : int
        Number of bundle points.
    alpha : float
        Weight of the bundle points.
    iterations : int
        Number of iterations of the algorithm.
    distance_metric : str
        Distance metric to use. Currently only 'euclidean' and 'mahalanobis' are supported.

    Returns
    -------
    points : array-like
        The final point distribution.
    """

    bounds = (0, 0.1, 1, 0.9)
    points = np.random.rand(num_points, 2)

    # x = np.linspace(0, 1, int(np.sqrt(num_points)))
    # y = np.linspace(0, 1, int(np.sqrt(num_points)))
    # X, Y = np.meshgrid(x, y)
    # points = np.column_stack([X.ravel(), Y.ravel()])

    bundle = 0.5 * np.random.rand(num_bundle, 2) + 0.25  # static for now

    for i in range(iterations):
        points = step_lloyd(
            points,
            bounds,
            bundle,
            alpha=alpha,
            distance_metric=distance_metric,
            beta_outside=1e16,
        )

    return points


def create_sphere(data):
    """
    Prend un tableau de points 2D et le transforme en un disque de centre (0.5, 0.5) et de rayon 0.4.
    Les points sont ensuite déplacés pour être centrés sur (0, 0) et dilatés pour avoir un rayon de 1.
    Les points en dehors de la sphère sont éliminés.

    Parameters
    ----------
    data : array-like
        Tableau de points 2D

    Returns
    -------
    data : array-like
        Tableau de points 2D transformés
    """

    centre_x = 0.5
    centre_y = 0.5
    rayon = 0.4
    distances = np.sqrt((data[:, 0] - centre_x) ** 2 + (data[:, 1] - centre_y) ** 2)
    data = data[distances <= rayon]
    data[:, 1] = data[:, 1] + 1
    data[:, 0] = 2 * data[:, 0]
    data[:, 1] = 2 * data[:, 1]

    return data


def filtrer_points_par_indicateur(array):
    """
    Filtrer les points d'un tableau en gardant uniquement ceux qui apparaissent
    au moins 4 fois dans la dernière colonne.

    Parameters
    ----------
    array : array-like
        Tableau de points 2D avec un indicateur en dernière colonne

    Returns
    -------
    array : array-like
        Tableau de points 2D filtrés
    """

    unique, counts = np.unique(array[:, -1], return_counts=True)
    occur_dict = dict(zip(unique, counts))

    array = np.array([p for p in array if occur_dict[p[-1]] >= 4])

    return array


def get_inset_voronoi_vertices(points, bounds, bw=0.1):
    """
    Compute the inset Voronoi vertices of a set of points within a given box.

    Parameters
    ----------
    points : array-like
        Input array of points
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box
    bw : float, optional
        Width of the inset region

    Returns
    -------
    inset_vertices : list of array-like
        List of arrays of coordinates of the inset Voronoi vertices
    num_cells : list of int
        List of indices of the regions corresponding to the inset vertices
    """

    xmin, xmax, ymin, ymax = bounds
    vor = Voronoi(points)
    inset_vertices = []
    num_cells = []

    for idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue
        poly_coords = np.array([vor.vertices[i] for i in region])
        if not np.all(
            (xmin <= poly_coords[:, 0])
            & (poly_coords[:, 0] <= xmax)
            & (ymin <= poly_coords[:, 1])
            & (poly_coords[:, 1] <= ymax)
        ):
            continue
        center = np.mean(poly_coords, axis=0)
        vecs = poly_coords - center
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        directions = vecs / norms
        adjusted = poly_coords - 0.5 * bw * directions
        inset_vertices.append(adjusted)
        num_cells.append(idx)

    return inset_vertices, num_cells


def get_inset_voronoi_vertices_option_boutplat(points, bounds, bw=0.1):
    """
    Compute the inset Voronoi vertices of a set of points within a given box.

    Parameters
    ----------
    points : array-like
        Input array of points
    bounds : tuple
        (min_x, min_y, max_x, max_y) coordinates of the bounding box
    bw : float, optional
        Width of the inset region

    Returns
    -------
    inset_vertices : list of array-like
        List of arrays of coordinates of the inset Voronoi vertices
    num_cells : list of int
        List of indices of the regions corresponding to the inset vertices
    """

    xmin, xmax, ymin, ymax = bounds
    ymin, ymax = 0.1, 0.9
    vor = Voronoi(points)
    box_poly = box(xmin, ymin, xmax, ymax)

    inset_vertices = []
    num_cells = []

    for idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue

        poly = Polygon([vor.vertices[i] for i in region])
        poly_clipped = poly.intersection(box_poly)
        if poly_clipped.is_empty or poly_clipped.geom_type != "Polygon":
            continue

        # Y : La on s'emmerde avec des buffers de connards
        poly_inset = poly_clipped.buffer(-0.5 * bw, join_style=2)
        if poly_inset.is_empty or poly_inset.geom_type != "Polygon":
            continue

        coords = list(poly_inset.exterior.coords)[:-1]
        inset_vertices.append(np.array(coords))
        num_cells.append(idx)

    return inset_vertices, num_cells


# run_lloyd()
# run_lloyd(distance_metric='bray_curtis')
# run_lloyd(distance_metric='canberra')
# run_lloyd(distance_metric='inf')
# run_lloyd(distance_metric='x_axis')
# run_lloyd(distance_metric='directional_distance')
# run_lloyd(distance_metric='y_axis')
# run_lloyd(distance_metric='mahalanobis')
# run_lloyd(distance_metric=0)
# run_lloyd(distance_metric=1)


def load_config(config_path):
    """
    Load a configuration from a JSON file.

    Parameters
    ----------
    config_path : str
        Path to the configuration file

    Returns
    -------
    config : dict
        Configuration dictionary
    """

    defaults = {
        "n_points": 1000,
        "n_bundle": 50,
        "alpha_hamiltonien": 0.2,
        "n_iter": 20,
        "dist_metric": "euclidean",
        "seed_number": 1234,
        "barWidth": 0.001,
        "bounds": [-0.5, 1.5, 0.0, 1.0],
    }

    with open(config_path, "r") as f:
        config = json.load(f)

    # Update only missing keys with defaults
    for key, value in defaults.items():
        if key not in config:
            config[key] = value

    return config


if __name__ == "__main__":
    # Load configuration from JSON file
    config = load_config(sys.argv[1])

    n_points = config["n_points"]
    seed_number = config["seed_number"]
    alpha_hamiltonien = config["alpha_hamiltonien"]
    barWidth = config["barWidth"]
    bounds = tuple(config["bounds"])  # (xmin, xmax, ymin, ymax)
    n_bundle = config["n_bundle"]
    n_iter = config["n_iter"]
    dist_metric = config["dist_metric"]

    ################################################

    np.random.seed(seed_number)

    points = run_lloyd(
        num_points=n_points,
        num_bundle=n_bundle,
        alpha=alpha_hamiltonien,
        iterations=n_iter,
        distance_metric=dist_metric,
    )

    # FIXME: il faudrait utiliser la version 'boutplat' que si un flag est activé
    # lhyphen_point, num_cells = get_inset_voronoi_vertices(points, bounds, bw=barWidth)
    lhyphen_point, num_cells = get_inset_voronoi_vertices_option_boutplat(
        points, bounds, bw=barWidth
    )

    lhyphen_format_x = []
    lhyphen_format_y = []
    lhyphen_format_cell = []
    for i in range(len(lhyphen_point)):
        masse_ponctuelle = lhyphen_point[i]

        for j in range(len(masse_ponctuelle)):
            if masse_ponctuelle[j][0] < 0.5 * (1 / 3) or masse_ponctuelle[j][0] > 1 - (
                0.5 * (1 / 3)
            ):
                continue
            # if masse_ponctuelle[j][1] < 0 or masse_ponctuelle[j][1] > 1:
            #    continue
            lhyphen_format_x.append(masse_ponctuelle[j][0])
            lhyphen_format_y.append(masse_ponctuelle[j][1])
            lhyphen_format_cell.append(num_cells[i])

    lhyphen_format_x = np.array(lhyphen_format_x)
    lhyphen_format_y = np.array(lhyphen_format_y)
    lhyphen_format_cell = np.array(lhyphen_format_cell)
    result = np.column_stack((lhyphen_format_x, lhyphen_format_y, lhyphen_format_cell))

    cell_ids = result[:, 2]
    unique, counts = np.unique(cell_ids, return_counts=True)
    occurrence_dict = dict(zip(unique, counts))
    valid_ids = {key for key, value in occurrence_dict.items() if value > 3}
    # valid_ids = {key for key, value in occurrence_dict.items() if value > 5}
    result = result[np.isin(cell_ids, list(valid_ids))]

    # result = create_sphere(result)
    # result = filtrer_points_par_indicateur(result)

    # tol = 1.01*min(result[:,1])
    # for i in range(len(result)):
    #    if result[i,1] > 1-tol:
    #        result[i,1] = result[i,1] + 0.02
    #    if result[i,1] < tol:
    #        result[i,1] = result[i,1] - 0.02

    np.savetxt("nodefile.txt", result, fmt="%1.16g %1.16g %d")
