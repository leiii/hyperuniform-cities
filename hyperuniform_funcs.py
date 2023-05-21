# -*- coding: utf-8 -*-
# Hyperuniformity
# author: Lei Dong
# contact: arch.dongl@gmail.com

import numpy as np
import statsmodels.api as sm

from scipy import stats
from scipy.spatial import distance, KDTree




def points_density_in_disk(pts, num_sample=10000, lmin=1., lmax=100., dl=1.):
    """Compute the point density vs. variation in 2D.
    
    Arguments:
        pts           points array
        lmin, lmax    min and max of the disk radius
        dl            length that the disk radius increases at each step
        num_sample    number of sampled disks under each radius

    Returns:
        sigma2        the list [radius, point density variation]
    """
    pts = np.array(pts)
    xmax, ymax = np.max(pts, axis=0)
    xmin, ymin = np.min(pts, axis=0)

    
    rhos = {}
    for l in np.arange(lmin, lmax + 1.1 * dl, dl):
        rhos[l] = []

        # Random sampling in the plane of (L-2l) * (L-2l)
        cxs = np.random.uniform(xmin+l, xmax-l, num_sample)
        cys = np.random.uniform(ymin+l, ymax-l, num_sample)
        
        for i in range(num_sample):
            X = np.array([[cxs[i], cys[i]]])
            dist_array = distance.cdist(X, pts, 'euclidean')
            density = (dist_array < l).sum() / (np.pi * l * l)
            rhos[l].append(density)
        
    # Calculate variation of density
    sigma2 = []
    for l in rhos:
        sigma2.append([l, np.var(rhos[l])])
    sigma2.sort()
    return sigma2




def structure_factor(pts, num_ks = 10000, num_bins = 100, scale = False):
    """Compute the structure factor as a function of wavevector k.
    We use the equation (2) given by Kai Zhang (arXiv:1606.03610v2):
    $S(k) = 1/N * <|\sum cos(k \dot r_i)|^2 + |\sum sin(k \dot r_i)|^2>$
    
    Arguments:
        pts           points array
        num_ks        the number of wavevectors k
        num_bins      the number of bins for the modulus k
        scale         scale the wavevectors to the system size

    Returns:
        sf_k          [the modulus k = |k_vec|, structure factor]
        sf_k_bins     [binned modulus k, average structure factor]
    """
    pts = np.array(pts)
    N = len(pts)
    xmax, ymax = np.max(pts, axis=0)
    xmin, ymin = np.min(pts, axis=0)

    L = min(xmax-xmin, ymax-ymin) # The size of the system
    dk = 2 * np.pi / L # delta wavelength

    # Two ways to generate wave vector
    size = int(np.sqrt(num_ks))

    if scale == False:
        scalefactor = 1
    else:
        scalefactor = L/size

    ks = {}
    for kx in range(1, size):
        for ky in range(1, size):
            k_vec = dk * np.array([kx * scalefactor, ky * scalefactor])
            sf_sin = 0
            sf_cos = 0
            for pt in pts:
                sf_sin += np.sin(k_vec.dot(pt))
                sf_cos += np.cos(k_vec.dot(pt))
            sf = ((sf_sin ** 2) + (sf_cos ** 2)) / N
            k_modulus = np.linalg.norm(k_vec)
            if k_modulus not in ks:
                ks[k_modulus] = [sf]
            else:
                ks[k_modulus].append(sf)

    sf_k = []
    for k in ks:
        sf_k.append([k, np.mean(ks[k])])
    sf_k.sort()


    # Aggregate data points into bins
    x = [i[0] for i in sf_k]
    y = [i[1] for i in sf_k]

    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y,
        statistic='mean', bins=num_bins)
    
    sf_k_bins = []
    for i in range(len(bin_means)):
        sf_k_bins.append([(bin_edges[i] + bin_edges[i+1])/2, bin_means[i]])


    return sf_k, sf_k_bins



def pair_correlation(pts, rmin=0.0, rmax=10.0, dr=0.1):
    """Compute the two-dimensional pair correlation function (i.e., the 
    radial distribution function).
    ref: https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation
    
    Arguments:
        pts                 points array
        rmin, rmax          min and max of the radius of disk
        dr                  increment for increasing radius of disk

    Returns: 
        g_average           an array containing the correlation function g(r)
        radii               an array containing the radii of the  disk
        interior_indices    indices of reference particles
    """
    pts = np.array(pts)
    xmax, ymax = np.max(pts, axis=0)
    xmin, ymin = np.min(pts, axis=0)
    density = len(pts) * 1.0 / (xmax-xmin) / (ymax-ymin)
    #print('The density of the systems is %s' % density)
    
    bools1 = pts > np.array([xmin+rmax, ymin+rmax])
    bools1 = np.all(bools1, axis=1)
    bools2 = pts < np.array([xmax-rmax, ymax-rmax])
    bools2 = np.all(bools2, axis=1)
    
    interior_indices, = np.where(bools1 * bools2)
    num_interior_pts = len(interior_indices)
    if num_interior_pts < 1:
        raise RuntimeError ("No particles found for which a circle of radius rmax\
                will lie entirely within a square of side length.  Decrease rmax\
                or increase the size of the square.")
        
    edges = np.arange(rmin, rmax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = np.zeros([num_interior_pts, num_increments])
    radii = np.zeros(num_increments)

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_pts):
        index = interior_indices[p]
        pt_select = np.array([pts[index]])
        d = np.sqrt(np.sum((pt_select - pts)**2, axis=1))
        
        d[index] = 2.0 * rmax # otherwise, the distance for this point = 0

        (result, bins) = np.histogram(d, bins=edges)
        g[p, :] = result/density

    # Average g(r) for all interior particles and compute radii
    g_average = np.zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = np.mean(g[:, i]) / (np.pi * (rOuter**2 - rInner**2))

    return radii, g_average



