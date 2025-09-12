from numpy import append,any,array,diff,interp,invert,quantile, \
                    sort,searchsorted,nansum,nanmax
from numpy.random import choice,normal
from scipy.special import gamma, gammainc, ndtr, ndtri

# Kaplan-Meier Curve
def km_curve(x,c):
    '''
    Name:
        km_curve

    Purpose:
        Compute the empirical non-parametric Kaplan-Meier survival function for
        a given data set, accounting for left censoring (lower limits). If data
        are right censored (upper limits), convert to left censoring before use.

    Arguments:
        :x (*np.ndarray*): 1xN data set for which to compute the Kaplan-Meier
                survival function
        :c (*np.ndarray*): 1xN array of integers or boolians indicating whether
                the corresponding element of `x` is censored (c=1 or c=True) or
                uncensored (c=0 or c=False)

    Returns:
        :km_x (*np.ndarray*): sorted array of `x` for the Kaplan-Meier curve
        :km_y (*np.ndarray*): Kaplan-Meier survival curve for `x`
    '''
    # check if censor indicators are booleans
    if c.dtype != bool:
        c = c==1
    u = invert(c)
    # sort uncensored values and add a terminating value at the end
    km_x = append(sort(x[u]),[max(x)//1+1])
    # number of cases
    n = len(x)
    # survival curve = 1 - number of cases to be lost / number of cases
    km_y = array( list( map( lambda i: \
      1-(len(x[u][(x[u]<km_x[i])])+len(x[c][(x[c]<km_x[i-1])]))/n,\
      range(1,len(km_x)) ) ) )
    # at the end, no cases remain except the upper limits,
    # which all eventually go away at t -> infty
    km_y[-1] = 0
    return km_x,append([1],km_y)

def km_eval(x0,x,c,x0_err=None):
    '''
    Name:
        km_eval

    Purpose:
        Evaluate the empirical non-parametric Kaplan-Meier survival function for
        a given data `x` set at some specified value `x0`. This evaluation
        provides the quantile(s) associated with `x0`, which are a statistical
        assessment of the null hypothesis that `x0` occurs within `x`. If the
        returned value `p_x` is outside the range of [0.01,0.99], the null
        hypothesis can be rejected with confidence. The Kaplan-Meier statistic
        is unique in that it computes a survival function which accounts for
        left censoring (lower limits). If data are right censored
        (upper limits), the user should convert to left censoring before use.

    Arguments:
        :x0 (*np.ndarray*): value(s) at which to evaluate the Kaplan-Meier curve
        :x (*np.ndarray*): 1xN data set for which to compute the Kaplan-Meier
                survival function
        :c (*np.ndarray*): 1xN array of integers or boolians indicating whether
                the corresponding element of `x` is censored (c=1 or c=True) or
                uncensored (c=0 or c=False)

    Keyword Arguments:
        :x0_err (*np.ndarray*): 1xN or 2xN array of uncertainties in `x0`. If
                provided, the Kaplan-Meier curve will also be evaluated at
                the uncertainties in `x0`. Default is `None`.

    Returns:
        :p_x (*np.ndarray*): Kaplan-Meier survival curve evaluated at `x0`,
                indicating the probability that `x0` is associated with `x`.
                Values of `p_x` < 0.01 or > 0.99 suggest rejection of the null
                hypthosesis that `x0` occurs within `x`.
        :p_e (*np.ndarray*): (optional) if `x0_err` is provided, a 1x2 array of
                uncertainties in `p_x` based on uncertainties in `x0`.
                Kaplan-Meier survival curve for `x` is evaluated at
                `x0`-`x0_err` and `x0`+`x0_err`, and the difference is returned.
    '''
    if any(x0_err==None):
        return interp(x0,*km_curve(x,c))
    if x0_err.ndim == 2:
        q = interp([x0-x0_err[0],x0,x0+x0_err[1]],*km_curve(x,c))
        return q[1],diff(q[::-1],axis=0).flatten()
    else:
        q = interp([x0-x0_err,x0,x0+x0_err],*km_curve(x,c))
        return q[1],diff(q[::-1],axis=0).flatten()

def km_var(x0,x,c,n_samp=1000,method='boot',xerr=None):
    '''
    Name:
        km_var

    Purpose:
        Determine the variations in the empirical non-parametric Kaplan-Meier
        survival function for a given data set `x` at a particular value `x0`,
        accounting for left censoring (lower limits). Variations are computed
        by either Monte Carlo simulation accounting for uncertainties `xerr` in
        `x` or by bootstrapping accounting for outliers in `x`. If data
        are right censored (upper limits), convert to left censoring before use.

    Arguments:
        :x0 (*float*): value at which to evaluate the Kaplan-Meier curve
        :x (*np.ndarray*): 1xN data set for which to compute the Kaplan-Meier
                survival function
        :c (*np.ndarray*): 1xN array of integers or boolians indicating whether
                the corresponding element of `x` is censored (c=1 or c=True) or
                uncensored (c=0 or c=False)

    Keyword Arguments:
        :n_samp (*int*): number of trials to compute, recommended no less than
                100. Default is 1000.
        :xerr (*np.ndarray*): 1xN array of uncertainties in `x`. Required if the
                `method` keyword is set to 'mc' or 'monte carlo'
        :method (*str*): string indicating the method to use for evalutating
                variations in the Kaplan-Meier curve. Options are Monte Carlo
                simulation ('mc' or 'monte carlo') to account for uncertainties
                in `x` or bootstrapping ('boot' or 'bootstrap') to account for
                outliers in `x`. Monte Carlo simulation requires passing 1xN
                array for keyword `xerr`. Default is 'boot'.

    Returns:
        :km_med (*float*): median Kaplan-Meier curve value for `x0`
        :km_sig (*np.ndarray*): 16-84 uncertainties in the Kaplan-Meier curve
                evaluated at `x0`.
    '''
    if 'boot' in method:
        ind = np.random.choice(len(x),(n_samp,len(x)))
        k = array(list(map(lambda i: interp(x0,*km_curve(x[i],c)),ind)))
    elif 'mc' in method or 'monte' in method:
        if any(xerr==None):
            print('Oops! Looks like you forgot to include your uncertainties!')
        X = normal(loc=x,scale=xerr,size=(n_samp,len(x)))
        k = array(list(map(lambda xi: interp(x0,*km_curve(xi,c)),X)))
    q = quantile(k,[0.1587,0.5,0.8413])
    return q[1],*diff(q)

# Mantel Log Rank Test for Kaplan Meier survival curves
def km_logrank(x1,c1,x2,c2,estimator='pike'):
    '''
    Name:
        km_logrank

    Purpose:
        Compute the log-rank test using Kaplan-Meier survival functions for
        two data sets, accounting for left censoring (lower limits). If data
        are right censored (upper limits), convert to left censoring before use.

    Arguments:
        :x1 (*np.ndarray*): 1xN data set for which to compute the Kaplan-Meier
                survival function; this is the reference data set
        :c1 (*np.ndarray*): 1xN array of integers or boolians indicating whether
                the corresponding element of `x1` is censored (c1=1 or c1=True) or
                uncensored (c1=0 or c1=False)
        :x2 (*np.ndarray*): 1xM data set for which to compute the Kaplan-Meier
                survival function; this is the test data set
        :c2 (*np.ndarray*): 1xM array of integers or boolians indicating whether
                the corresponding element of `x2` is censored (c2=1 or c2=True) or
                uncensored (c2=0 or c2=False)

    Keyword Arguments:
        :estimator (*str*): method for computing the estimator of the variance-
                weighted square difference between distributions. Options are
                `mantle` for Mantel 1966, or `pike` for Pike 1972.
                Default is `pike`, the more conservative (Peto & Pike 1973).

    Returns:
        :D (*float*): maximum difference between the two survival curves
        :T (*float*): variance-weighted squared difference between the two
                survival curves
        :p (*float*): probability that the null hypothesis is true assuming that
                the variance-weighted squared difference is drawn from a chi_sq
                distribution with one degree of freedom
    '''
    # intervals
    ivls  = sort(append(x1,x2))
    # total events, censored and uncensored, up to each point in full sample
    Y1    = searchsorted(sort(x1),ivls).astype(float)
    Y2    = searchsorted(sort(x2),ivls).astype(float)
    # events which are uncensored
    u1    = c1==False
    u2    = c2==False
    # total uncensored events up to each point in full sample
    d1    = searchsorted(sort(x1[u1]),ivls).astype(float)
    d2    = searchsorted(sort(x2[u2]),ivls).astype(float)
    # sum up for both data sets
    Y     = Y1+Y2
    d     = d1+d2
    # expected uncensored events in new/test data set
    # based on the original/reference data set
    E1    = d*Y1/Y
    E2    = d*Y2/Y
    # variance based on the two data sets
    V     = Y1*Y2*d*(Y-d)/(Y**2.*(Y-1.))
    # variance-weighted differences
    if estimator.lower() not in ['mantel','pike']:
        raise ValueError(f'Estimator \'{estimator}\' not recognized.\n'+\
                'Options are \'mantel\' or \'pike\'.')
    # Mantel 1966 estimator for T
    elif estimator.lower() == 'mantel' :
        T    = nansum(d1-E1)**2/nansum(V)
    # Pike 1972 estimator for T
    elif estimator.lower() == 'pike' :
        T    = nansum(d1-E1)**2/nansum(E1)+nansum(d2-E2)**2/nansum(E2)
    # p-value with Z^2 drawn from chi^2(1)
    pval  = 1.-gammainc(0.5,0.5*T)
    return abs(nanmax(d1-E1)),T,pval
