# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import glob


def ecdf(data):
    """
    Computes the empirical cumulative distribution function for a collection of provided data.

    Parameters
    ----------
    data : 1d-array, Pandas Series, or list
        One-dimensional collection of data for which the ECDF will
        be computed

    Returns
    -------
    x, y : 1d-arrays
        The sorted x data and the computed ECDF
    """
    return np.sort(data), np.arange(0, len(data)) / len(data)


def _log_prior_trace(trace, model):
    """
    Computes the contribution of the log prior to the log posterior.

    Parameters
    ----------
    trace : PyMC3 trace object.
        Trace from the PyMC3 sampling.
    model : PyMC3 model object
        Model under which the sampling was performed

    Returns
    -------
    log_prior_vals : nd-array
        Array of log-prior values computed elementwise for each point in the
        trace.

    Notes
    -----
    This function was modified from one produced by Justin Bois.
    http://bebi103.caltech.edu
    """
    # Iterate through each trace.
    try:
        points = trace.points()
    except:
        points = trace

    # Get the unobserved variables.
    priors = [var.logp for var in model.unobserved_RVs if type(var) == pm.model.FreeRV]

    def logp_vals(pt):
        if len(model.unobserved_RVs) == 0:
            return pm.theanof.floatX(np.array([]), dtype="d")

        return np.array([logp(pt) for logp in priors])

    # Compute the logp for each value of the prior.
    log_prior = (logp_vals(pt) for pt in points)
    return np.stack(log_prior)


def _log_post_trace(trace, model):
    r"""
    Computes the log posterior of a PyMC3 sampling trace.

    Parameters
    ----------
    trace : PyMC3 trace object
        Trace from MCMC sampling
    model: PyMC3 model object
        Model under which the sampling was performed.

    Returns
    -------
    log_post : nd-array
        Array of log posterior values computed elementwise for each point in
        the trace

    Notes
    -----
    This function was modified from one produced by Justin Bois
    http://bebi103.caltech.edu
    """

    # Compute the log likelihood. Note this is improperly named in PyMC3.
    log_like = pm.stats._log_post_trace(trace, model).sum(axis=1)

    # Compute the log prior
    log_prior = _log_prior_trace(trace, model)

    return log_prior.sum(axis=1) + log_like


def trace_to_dataframe(trace, model):
    r"""
    Converts a PyMC3 sampling trace object to a pandas DataFrame

    Parameters
    ----------
    trace, model: PyMC3 sampling objects.
        The MCMC sampling trace and the model context.

    Returns
    -------
    df : pandas DataFrame
        A tidy data frame containing the sampling trace for each variable  and
        the computed log posterior at each point.
    """

    # Use the Pymc3 utilitity.
    df = pm.trace_to_dataframe(trace)

    # Include the log prop
    df["logp"] = _log_post_trace(trace, model)
    return df


def compute_statistics(df, varnames=None, logprob_name="logp"):
    r"""
    Computes the mode, hpd_min, and hpd_max from a pandas DataFrame. The value
    of the log posterior must be included in the DataFrame.
    """

    # Get the vars we care about.
    if varnames is None:
        varnames = [v for v in df.keys() if v is not "logp"]

    # Find the max of the log posterior.
    ind = np.argmax(df[logprob_name].values)

    # Instantiate the dataframe for the parameters.
    stat_df = pd.DataFrame(
        [], columns=["parameter", "mean", "median", "mode", "hpd_min", "hpd_max"]
    )
    for v in varnames:
        mode = df.iloc[ind][v]
        median = df[v].median()
        mean = df[v].mean()
        hpd_min, hpd_max = compute_hpd(df[v].values, mass_frac=0.95)
        stat_dict = dict(
            parameter=v,
            median=median,
            mean=mean,
            mode=mode,
            hpd_min=hpd_min,
            hpd_max=hpd_max,
        )
        stat_df = stat_df.append(stat_dict, ignore_index=True)

    return stat_df

def bin_by_events(df, bin_size, sortby='summed', average=['summed', 'fluct']):
    """
    Bins a given data set by number of events rather than bin width and
    computes the mean value of desired parameters.

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe containing data to bin and average.
    bin_size : int
        Number of events to consider for one bin.
    sortby : str
        The name of the column to sort the values by. Default is 'summed'
    average : list
        The quantities over which to average. These will be returned in the
        order they are provided.

    Returns
    -------
    average_vals : list with shape of `average`.
        The average of the quantities in each bin. This is a list of lists.
    """
    num_quantities = len(average)

    # Sort the dataframe.
    sorted_df = df.sort_values(sortby)

    # Set the bins.
    bins = np.arange(0, len(sorted_df) + bin_size, bin_size)
    averages = {i: np.zeros(len(bins)) for _, i in enumerate(average)}
    for k in average:
        averages[f'{k}_sem'] = np.zeros(len(bins))
    # Iterate through each bin and compute the average quanities.
    for i in range(1, len(bins)):
        # Slice the data frame.
        d = sorted_df.iloc[bins[i - 1]:bins[i]][average]
        for k in d.keys():
            val = np.mean(d[k].values)
            averages[f'{k}_sem'][i-1] = np.std(d[k].values) / np.sqrt(len(d))
            averages[k][i-1] = val
    return averages


def compute_hpd(trace, mass_frac):
    r"""
    Returns highest probability density region given by
    a set of samples.

    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For hreple, `massfrac` = 0.95 gives a
        95% HPD.

    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD

    Notes
    -----
    We thank Justin Bois (BBE, Caltech) for developing this function.
    http://bebi103.caltech.edu/2015/tutorials/l06_credible_regions.html
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)

    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)

    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[: n - n_samples]

    # Pick out minimal interval
    min_int = np.argmin(int_width)

    # Return interval
    return np.array([d[min_int], d[min_int + n_samples]])


def compute_mean_sem(df):
    """
    Computes the mean and standard error of the fold-change given a
    grouped pandas Series.
    """
    # Compute the properties
    mean_fc = df["fold_change"].mean()
    sem_fc = df["fold_change"].std() / np.sqrt(len(df))

    # Assemble the new pandas series and return.
    samp_dict = {"mean": mean_fc, "sem": sem_fc}
    return pd.Series(samp_dict)
