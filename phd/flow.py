import numpy as np
import fcsparser
import pandas as pd
from ._fit_bivariate_normal_AstroML import fit_bivariate_normal
import scipy.stats


# #######################
# Automated Gating
# #######################
def fit_2D_gaussian(df, x_val="FSC-H", y_val="SSC-H", log=False):
    """
    This function hacks astroML fit_bivariate_normal to return the mean
    and covariance matrix when fitting a 2D gaussian fuction to the data
    contained in the x_val and y_val columns of the DataFrame df.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not

    Returns
    -------
    mu : tuple.
        (x, y) location of the best-fit bivariate normal
    cov : 2 x 2 array
        covariance matrix.
        cov[0, 0] = variance of the x_val column
        cov[1, 1] = variance of the y_val column
        cov[0, 1] = cov[1, 0] = covariance of the data
    """
    if log:
        x = np.log10(df[x_val])
        y = np.log10(df[y_val])
    else:
        x = df[x_val]
        y = df[y_val]

    # Fit the 2D Gaussian distribution using atroML function
    mu, sigma_1, sigma_2, alpha = fit_bivariate_normal(x, y, robust=True)

    # compute covariance matrix from the standar deviations and the angle
    # that the fit_bivariate_normal function returns
    sigma_xx = (sigma_1 * np.cos(alpha)) ** 2 + (sigma_2 * np.sin(alpha)) ** 2
    sigma_yy = (sigma_1 * np.sin(alpha)) ** 2 + (sigma_2 * np.cos(alpha)) ** 2
    sigma_xy = (sigma_1 ** 2 - sigma_2 ** 2) * np.sin(alpha) * np.cos(alpha)

    # put elements of the covariance matrix into an actual matrix
    cov = np.array([[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]])

    return mu, cov


# #################
def gauss_interval(df, mu, cov, x_val="FSC-H", y_val="SSC-H", log=False):
    """
    Computes the of the statistic

    (x - µx)'Σ(x - µx)

    for each of the elements in df columns x_val and y_val.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    mu : array-like.
        (x, y) location of bivariate normal
    cov : 2 x 2 array
        covariance matrix
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not.

    Returns
    -------
    statistic_gauss : array-like.
        array containing the result of the linear algebra operation:
        (x - µx)'sum(x - µx)
    """
    # Determine that the covariance matrix is not singular
    det = np.linalg.det(cov)
    if det == 0:
        raise NameError("The covariance matrix can't be singular")

    # Compute the vector x defined as [[x - mu_x], [y - mu_y]]
    if log is True:
        x_vect = np.log10(np.array(df[[x_val, y_val]]))
    else:
        x_vect = np.array(df[[x_val, y_val]])

    x_vect[:, 0] = x_vect[:, 0] - mu[0]
    x_vect[:, 1] = x_vect[:, 1] - mu[1]

    # compute the inverse of the covariance matrix
    inv_sigma = np.linalg.inv(cov)

    # compute the operation
    interval_array = np.zeros(len(df))
    for i, x in enumerate(x_vect):
        interval_array[i] = np.dot(np.dot(x, inv_sigma), x.T)

    return interval_array


def gaussian_gate(df, alpha, x_val="FSC-A", y_val="SSC-A", log=True, verbose=False):
    """
    Function that applies an "unsupervised bivariate Gaussian gate" to the data
    over the channels x_val and y_val.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    alpha : float. [0, 1]
        fraction of data aimed to keep. Used to compute the chi^2 quantile
        function
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not
    verbose : bool.
        indicate if the percentage of data kept should be print

    Returns
    -------
    df_thresh : DataFrame
        Pandas data frame to which the automatic gate was applied.
    """

    # Perform sanity checks.
    if alpha < 0 or alpha > 1:
        return RuntimeError("`alpha` must be a float between 0 and 1.")

    data = df[[x_val, y_val]]
    # Fit the bivariate Gaussian distribution
    mu, cov = fit_2D_gaussian(data, log=log, x_val=x_val, y_val=y_val)

    # Compute the statistic for each of the pair of log scattering data
    interval_array = gauss_interval(data, mu, cov, log=log, x_val=x_val, y_val=y_val)

    # Find which data points fall inside the interval
    idx = interval_array <= scipy.stats.chi2.ppf(alpha, 2)

    # print the percentage of data kept
    if verbose:
        print(
            """
        with parameter alpha={0:0.2f}, percentage of data kept = {1:0.2f}
        """.format(
                alpha, np.sum(idx) / len(df)
            )
        )
    return df[idx]


# #######################
# File Parsing Utilities
# #######################


def fcs_to_csv(path, file_name, save_metadata=True):
    r"""
    Reads in a Flow Cytometry Standard (FCS) file and exports all content
    directly to an easily parseable csv fie.

    Parameters
    ----------
    path : str
        Path to .fcs file
    file_name : str
        Path to save file to .csv
    save_metadata : bool
        If True, a metadata file will also be saved. It will have the name of
        `path` with `_metadata.csv`
    """

    # Ensure provided file is actually .fcs
    if path.split(".")[-1] is not ".fcs":
        raise RuntimeError("`path` is not an FCS file.")

    meta, data = fcsparser.parse(path)
    data.to_csv(file_name, index=False)

    if save_metadata:
        meta_df = pd.DataFrame(meta)
        meta_name = "{0}_metadata.csv".format(path[:-4])
        meta_df.to_csv(meta_name, index=False)
