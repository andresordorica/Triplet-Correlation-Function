import numpy as np
cimport numpy as np
from scipy.spatial.distance import pdist, squareform
import itertools
from libc.math cimport acos, cos, sqrt
from tqdm import tqdm
from cython.view cimport array as cvarray



def _tcf_self(np.ndarray[np.float64_t, ndim=3] positions_array, double tolerance=0.010, double angle=60):
    cdef Py_ssize_t n_frames, N, f
    cdef double cos_, d_1, d_2, d_jk, factor, cos_angle, angle_deg
    cdef Py_ssize_t i, j, k
    cdef np.ndarray[np.float64_t, ndim=2] array
    cdef np.ndarray[np.float64_t, ndim=2] d
    cdef Py_ssize_t n_bins
    cdef np.ndarray[np.float64_t, ndim=1] valid_distances
    cdef np.ndarray[np.float64_t, ndim=2] valid_angles
    cdef np.ndarray[np.int64_t, ndim=1] g_
    cdef np.ndarray[np.float64_t, ndim=1] edges
    cdef list an_
    cdef np.ndarray[np.float64_t, ndim=1] an
    cdef int bin_width


    n_frames = positions_array.shape[0]
    N = positions_array.shape[1]
    cos_ = (1 - cos(np.radians(angle))) * 2
    valid_distances = np.empty(0, dtype=np.float64)
    bin_width  = 2
    r_range = np.array([20.0, 180.0], dtype=np.float64)
    n_bins = int((r_range[1] - r_range[0]) / bin_width)
    valid_angles = np.empty((n_frames, n_bins), dtype=np.float64)
    pbar_frames = tqdm(total=n_frames, desc='Frames')
   
    for f in range(n_frames):
        pbar_frames.update(1)
        array = positions_array[f]
        N = array.shape[0]
        d = squareform(pdist(array, metric='euclidean'))
        an_ = []
        combinations = list(itertools.combinations(range(N), 3))
        i_indices, j_indices, k_indices = zip(*combinations)

        valid_i_indices = []
        valid_j_indices = []
        valid_k_indices = []
        valid_d_1_elements = []
        pbar_inner = tqdm(total=len(combinations), desc='Combinations', leave=False)
        for idx in range(len(i_indices)):
            pbar_inner.update(1)
            i = i_indices[idx]
            j = j_indices[idx]
            k = k_indices[idx]
            d_1 = abs(d[i, j] - d[i, k])
            
            if d_1 < tolerance:
                if d_1 < 0.34:
                    cos_angle = (d[i, j] ** 2 + d[i, k] ** 2 - d[j, k] ** 2) / (2 * d[i, j] * d[i, k])
                    angle_deg = np.degrees(acos(cos_angle))
                    an_.append(angle_deg)

                factor = sqrt((d[i, j] ** 2) * cos_)
                d_2 = abs(factor - d[j, k])
                
                if abs(d_1 - d_2) < tolerance:
                    valid_i_indices.append(i)
                    valid_j_indices.append(j)
                    valid_k_indices.append(k)
                    valid_d_1_elements.append(d[i, j])

        valid_i_indices = np.array(valid_i_indices)
        valid_j_indices = np.array(valid_j_indices)
        valid_k_indices = np.array(valid_k_indices)
        valid_d_1_elements = np.array(valid_d_1_elements)
        pbar_inner.close()  # Close progress bar for frames loop
        an = np.array(an_, dtype=np.float64)
        valid_distances = np.concatenate((valid_distances, valid_d_1_elements))
        g_, edges = np.histogram(an, range=r_range, bins=n_bins)
        valid_angles[f] = g_

        

    pbar_frames.close()  
    return valid_distances, valid_angles




def _tcf_(np.ndarray[np.float64_t, ndim=3] positions_array,int A, double tolerance=0.010, double angle=60):
    cdef Py_ssize_t n_frames, N, f
    cdef double cos_, d_1, d_2, d_jk, factor, cos_angle, angle_deg
    cdef Py_ssize_t i, j, k
    cdef np.ndarray[np.float64_t, ndim=2] array
    cdef np.ndarray[np.float64_t, ndim=2] d
    cdef Py_ssize_t n_bins
    cdef np.ndarray[np.float64_t, ndim=1] valid_distances
    cdef np.ndarray[np.float64_t, ndim=2] valid_angles
    cdef np.ndarray[np.int64_t, ndim=1] g_
    cdef np.ndarray[np.float64_t, ndim=1] edges
    cdef list an_
    cdef np.ndarray[np.float64_t, ndim=1] an
    cdef int bin_width


    n_frames = positions_array.shape[0]
    N = positions_array.shape[1]
    cos_ = (1 - cos(np.radians(angle))) * 2
    valid_distances = np.empty(0, dtype=np.float64)
    bin_width  = 2
    r_range = np.array([20.0, 180.0], dtype=np.float64)
    n_bins = int((r_range[1] - r_range[0]) / bin_width)
    valid_angles = np.empty((n_frames, n_bins), dtype=np.float64)
    pbar_frames = tqdm(total=n_frames, desc='Frames')
  
    for f in range(n_frames):
        pbar_frames.update(1)
        array = positions_array[f]
        N = array.shape[0]
        d = squareform(pdist(array, metric='euclidean'))
        an_ = []
        combinations_ = list(itertools.combinations(range(A), 1))
        extracted_elements = [t[0] for t in combinations_]
        combinations2 = list(itertools.combinations(range(A, N), 2))
        combinations = [(element,) + t for element in extracted_elements for t in combinations2]
        i_indices, j_indices, k_indices = zip(*combinations)

        valid_i_indices = []
        valid_j_indices = []
        valid_k_indices = []
        valid_d_1_elements = []
        pbar_inner = tqdm(total=len(combinations), desc='Combinations', leave=False)
        for idx in range(len(i_indices)):
            pbar_inner.update(1)
            i = i_indices[idx]
            j = j_indices[idx]
            k = k_indices[idx]
            d_1 = abs(d[i, j] - d[i, k])
            
            if d_1 < tolerance:
                if d_1 < 0.34:
                    cos_angle = (d[i, j] ** 2 + d[i, k] ** 2 - d[j, k] ** 2) / (2 * d[i, j] * d[i, k])
                    angle_deg = np.degrees(acos(cos_angle))
                    an_.append(angle_deg)

                factor = sqrt((d[i, j] ** 2) * cos_)
                d_2 = abs(factor - d[j, k])
                
                if abs(d_1 - d_2) < tolerance:
                    valid_i_indices.append(i)
                    valid_j_indices.append(j)
                    valid_k_indices.append(k)
                    valid_d_1_elements.append(d[i, j])

        valid_i_indices = np.array(valid_i_indices)
        valid_j_indices = np.array(valid_j_indices)
        valid_k_indices = np.array(valid_k_indices)
        valid_d_1_elements = np.array(valid_d_1_elements)
        an = np.array(an_, dtype=np.float64)
        valid_distances = np.concatenate((valid_distances, valid_d_1_elements))
        g_, edges = np.histogram(an, range=r_range, bins=n_bins)
        pbar_inner.close()
        valid_angles[f] = g_
        

    pbar_frames.close()
    return valid_distances, valid_angles


