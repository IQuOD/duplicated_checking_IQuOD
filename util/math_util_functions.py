import numpy as np
import pandas as pd
import os
import scipy.io as sio
from scipy.stats import zscore

def process_file_pairs_generic(file_pairs):
    from itertools import combinations
    unique_pairs = set()  ###  unique the list

    for pair in file_pairs:
        # all possible double pairs (sepreate the triple pairs etc.)
        for combo in combinations(pair, 2):
            sorted_combo = tuple(sorted(combo))
            unique_pairs.add(sorted_combo)

    return list(unique_pairs)


def normalizeData(data, new_min, new_max):
    """
    Normalize data to a new range [new_min, new_max].

    Parameters:
    data : numpy.ndarray
        Input data.
    new_min : float
        New minimum value of the range.
    new_max : float
        New maximum value of the range.

    Returns:
    numpy.ndarray
        Normalized data.
    """
    # Calculate the minimum and maximum values of the data, ignoring NaNs
    min_data = np.nanmin(data, axis=0)
    max_data = np.nanmax(data, axis=0)

    # Initialize normalized array with the same shape as input data
    normalized = np.zeros_like(data, dtype=float)

    # Normalize data where the min and max values are not equal
    valid_scale = max_data != min_data
    normalized[:, valid_scale] = (data[:, valid_scale] - min_data[valid_scale]) / (max_data[valid_scale] - min_data[valid_scale])

    # Scale to [new_min, new_max]
    normalized[:, valid_scale] = normalized[:, valid_scale] * (new_max - new_min) + new_min

    # Where min and max are the same, set to new_min (or any specific rule)
    normalized[:, ~valid_scale] = new_min

    return normalized

def standardizeData(data, targetMean, targetStd):
    """
    Standardize the data to a specified mean and standard deviation.
    
    Parameters:
        data (numpy.ndarray): The input data.
        targetMean (float): The target mean.
        targetStd (float): The target standard deviation.
        
    Returns:
        numpy.ndarray: The standardized data adjusted to the target mean and standard deviation.
    """
    # Calculate the mean and standard deviation of the data, ignoring NaNs

    meanData = np.nanmean(data,0)

    stdData = np.nanstd(data,0,ddof=1)


    # Normalize data to mean 0 and std 1
    standardizedData = (data - meanData) / stdData
    
    # Scale and shift the normalized data to the target mean and std
    standardizedData = standardizedData * targetStd + targetMean
    
    return standardizedData


def entropy_weight(DNA_series):
    DNA_mapped = normalizeData(DNA_series, 0, 1)
    DNA_mapped[DNA_series == 0] = 0

    B = DNA_mapped.copy()
    B[B == 0] = 0.00001
    B[B == 1] = 0.99999

    n, m = B.shape  # n samples, m variables
    p = np.zeros_like(B)  # proportion of ith sample in the jth variable
    dd = np.nansum(B, axis=0)  # sum across columns, ignoring NaNs

    for j in range(m):
        p[:, j] = B[:, j] / dd[j]

    # Calculate the entropy of the jth variable e(j)
    k = 1 / np.log(n)
    e = np.zeros(m)
    for j in range(m):
        e[j] = -k * np.nansum(p[:, j] * np.log(p[:, j]))

    # Calculate information entropy redundancy
    d = 1 - e
    weight = d / np.nansum(d)  # calculate the weight

    # Calculate composite score
    s = 100 * np.dot(weight, B.T)
    return weight, s

def PCA_DNA_profiles(DNA_series,threshold=0.9):
    """
    Perform PCA on DNA profile data and reduce dimensionality based on 90% variance threshold.

    Parameters:
    DNA_series_copy (numpy.ndarray): The input DNA series data.

    Returns:
    numpy.ndarray: The reduced DNA profile data based on PCA.
    """
    import numpy as np
    from scipy.stats import zscore
    
    z = zscore(DNA_series,axis=0)
    z = np.nan_to_num(z, nan=0.0)

    # Compute covariance matrix
    R = np.cov(z, rowvar=False)

    # Calculate the eigenvectors and eigenvalues of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(R)

    # Sort the eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    # Determine the number of components to retain based on 90% variance
    cumulative_variance = np.cumsum(sorted_eigenvalues)
    total_variance = cumulative_variance[-1]
    num_components = np.where(cumulative_variance / total_variance >= threshold)[0][0] + 1

    # Project the normalized data onto the retained eigenvectors
    DNA_new = np.dot(z, sorted_eigenvectors[:, :num_components])

    return DNA_new