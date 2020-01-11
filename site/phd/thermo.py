"""
A module for computing properties of various transcriptional
regulatory architectures.
"""
import numpy as np
import scipy.optimize


class MWC(object):
    r"""
    A base class for the Monod - Wyman - Changeux model for
    allostery.
    """

    def __init__(
        self,
        effector_conc=None,
        ka=None,
        ki=None,
        ep_ai=None,
        n_sites=2,
        log_transform=False,
    ):
        """
        Parameters
        ----------
        ep_ai : int, float, or array
            Difference in energy between the active and inactive allosteric
            states of the repressor. This should be in units of k_BT.
        ka, ki : ints, floats, or arrays
            The effector dissociation constants for the acitve and inactive
            state of the repressor.
        log_transform:  bool
            If True, the provided ka and ki are the log transform and will be
            exponentiated in the calculation of pact.
        effector_conc: int, float, or array
            Concentration of the allosteric effector molecule.
        n_sites : int, float or array
            Number of cooperative effector binding sites on the repressor.
            Default value is 2.
        """
        kwargs = dict(
            effector_conc=effector_conc, ka=ka, ki=ki, ep_ai=ep_ai, n_sites=n_sites
        )

        # Ensure values are provided.
        for k in kwargs.keys():
            if type(kwargs[k]) is None:
                raise RuntimeError("{0} is NoneType and must be defined.".format(k))

        # Assign the variables.
        self.c = effector_conc
        self.ep_ai = ep_ai
        self.n = n_sites
        if log_transform is True:
            self.ka = np.exp(ka)
            self.ki = np.exp(ki)
        else:
            self.ka = ka
            self.ki = ki

        # Ensure ka and ki are not zero.
        if type(ka) is float or int:
            _ka = np.array([ka])
        if type(ki) is float or int:
            _ki = np.array([ki])

        if (_ka == 0).any() or (_ki == 0).any():
            raise ValueError("ka and/or ki cannot be zero.")

        # Ensure positivity of values.
        positive_kwargs = dict(
            effector_conc=self.c, ka=self.ka, ki=self.ki, n_sites=self.n
        )
        for k in positive_kwargs.keys():
            val = positive_kwargs[k]
            if type(val) is float or int:
                val = np.array([val])
            if (val < 0).any():
                raise RuntimeError("{0} must be positive.".format(k))

    def pact(self):
        r"""
        Compute the probability of the active state at each provided parameter
        value

        Returns
        -------
        p_active : float or nd-array
            The probability of the active state evaluated at each value of
            effector_conc, ka, ki, and n_sites
        """
        c = self.c
        n = self.n
        ka = self.ka
        ki = self.ki
        ep_ai = self.ep_ai
        numer = (1 + c / ka) ** n
        denom = numer + np.exp(-ep_ai) * (1 + c / ki) ** n
        return numer / denom

    def saturation(self):
        r"""
        Computes the probability of the active state in the limit of
        saturating effector concentration.

        Returns
        -------
        saturation : float or nd-array
            Saturation value at each provided value of ka, ki, ep_ai, and
            n_sites.
        """
        ka = self.ka
        ki = self.ki
        ep_ai = self.ep_ai
        n = self.n
        return (1 + np.exp(-ep_ai) * (ka / ki) ** n) ** -1

    def leakiness(self):
        r"""
        COmputes the probability of the active state in the limit of zero effector.
        """
        return (1 + np.exp(-self.ep_ai)) ** -1


class SimpleRepression(object):
    r"""
    A base class for simple repression with an allosteric
    repressor.
    """

    def __init__(self, R, ep_r, n_ns=4.6e6, **kwargs):
        r"""
        Instantiate the SimpleRepression object.

        Parameters
        ----------
        R : int, float, or array
            Number of repressors in the system (per cell).
        ep_r : int, float or array
            Repressor-DNA binding energy in units of k_BT.
        n_ns : int or float
            Number of nonspecific DNA binding sites for the
            repressor molecule.
            Default value is the approximate length of the *E.
            coli* genome, 4.6e6 bp.
        **kwargs : dict or tuple
            kwargs for allosteric transcription factors see `MWC`
            documentation for more information.
        """
        # Define the variables.
        self.R = R
        self.ep_r = ep_r
        self.n_ns = n_ns

        # Ensure values make sense.
        positive_args = dict(R=R, n_ns=n_ns)
        for p in positive_args.keys():
            val = positive_args[p]
            if type(val) is float or int:
                val = np.array([val])
            if (val < 0).any():
                raise RuntimeError("{0} must be positive.".format(p))

        # Determine if transcription factor is allosteric
        if kwargs:
            self.allo = True
            self.mwc = MWC(**kwargs)
        else:
            self.allo = False

    def fold_change(self, wpa=True, num_pol=None, ep_pol=None, pact=False):
        r"""
        fold - change for simple repression.

        Parameters
        ----------
        wpa: bool
            If True, the weak promoter approximation is made and the state of
            polymerase being bound to the promoter is ignored.
        num_pol: int, float, or array
            Number of RNA Polymerase units per cell. This is required if
            `wpa == True`.
        ep_pol: int, float, or array
            RNAP - DNA binding energy in units of k_BT. This required if
            `wpa == True`.
        pact : float or array
            The probability of having an active repressor. If None is
            provided, the probability will be computed given effector_conc.

        Returns
        -------
        fold_change: float or nd - array
            Fold - change in gene expression evaluated at each value of c.
        """
        if self.allo == False:
            pact = 1
        else:
            if type(pact) == bool:
                pact = self.mwc.pact()
        # Compute repression and return inverse.
        repression = 1 + pact * (self.R / self.n_ns) * np.exp(-self.ep_r)
        return repression ** -1

    def saturation(self, wpa=True, num_pol=None, ep_pol=0):
        r"""
        Computes the fold - change in gene expression under saturating
        concentrations of effector. This function  is only defined for
        allosteric repressors.

        Parameters
        ----------
        wpa : bool
            If True, the weak promoter approximation will be applied.
        num_pol : int, float, or array
            The number of RNA Polymerase molecules per cell. This is required
            if `wpa == False`.
        ep_pol : int, float, or array
            The RNAP-DNA binding energy in units of k_BT. This is required if
            `wpa == False`
        Returns
        -------
        saturation: float or array
            The leakiness of the simple repression architecture.

        """
        if self.allo is False:
            raise RuntimeError(
                """Saturation is only defined for allosteric molecules. (`allosteric = True`)"""
            )
        # Compute the pact in limit of c -> inf.
        pact = self.mwc.saturation()
        return self.fold_change(wpa, num_pol, ep_pol, pact)

    def leakiness(self, wpa=True, num_pol=None, ep_pol=0):
        r"""
        Computes the fold-change in gene expression under a zero concentration
        of effector.

        Parameters
        ----------
        wpa : bool
            If True, the weak promoter approximation will be applied.
        num_pol : int, float, or array
            The number of RNA Polymerase molecules per cell. This is required
            if `wpa == False`.
        ep_pol : int, float, or array
            The RNAP-DNA binding energy in units of k_BT. This is required if
            `wpa == False`
        Returns
        -------
        leakiness: float or array
            The leakiness of the simple repression architecture.
        """
        # Compute the pact in the limit of c -> 0.
        if self.allo is True:
            pact = self.mwc.leakiness()
        else:
            pact = 1
        return self.fold_change(wpa, num_pol, ep_pol, pact)

    def dynamic_range(self, wpa=True, num_pol=None, ep_pol=0):
        r"""
        The dynamic range of the fold - change in response to an effector
        molecule. This property is only defined for allosteric molecules.

        Parameters
        ----------
        wpa : bool
            If True, the weak promoter approximation will be applied.
        num_pol : int, float, or array
            The number of RNA Polymerase molecules per cell. This is required
            if `wpa == False`.
        ep_pol : int, float, or array
            The RNAP-DNA binding energy in units of k_BT. This is required if
            `wpa == False`
        Returns
        -------
        dynamic_range: float or array
            The leakiness of the simple repression architecture.
        """
        # Compute the saturation and leakiness.
        sat = self.saturation(wpa, num_pol, ep_pol)
        leak = self.leakiness(wpa, num_pol, ep_pol)
        return sat - leak

    def ec50(self):
        """Computes the EC50 for allosteric architectures"""
        if self.allo is False:
            raise RuntimeError("EC50 defined only for allosteric architectures.")
        # Determine the user provided inputs.
        R = self.R
        n_ns = self.n_ns
        ep_r = self.ep_r
        ep_ai = self.mwc.ep_ai
        ka = self.mwc.ka
        ki = self.mwc.ki
        n_sites = self.mwc.n
        # Break it into pieces
        repression = 1 + (R / n_ns) * np.exp(-ep_r)
        numer = repression + (ka / ki) ** n_sites * (2 * np.exp(-ep_ai) + repression)
        denom = 2 * repression + np.exp(-ep_ai) + (ka / ki) ** n_sites * np.exp(-ep_ai)

        # Assemble the pieces of the ec50 calculation.
        ec50_numer = (ka / ki) - 1
        ec50_denom = (ka / ki) - (numer / denom) ** (1 / n_sites)
        return ka * ((ec50_numer / ec50_denom) - 1)

    def effective_hill(self):
        """Computes the effective hill coefficient of an allosteric repressor."""
        if self.allo == False:
            return RuntimeError(
                "Effective hill only defined for allosteric architectures"
            )

        # Define the parameters
        c = self.ec50()
        ka = self.mwc.ka
        ki = self.mwc.ki
        ep_ai = self.mwc.ep_ai
        n_sites = self.mwc.n
        R = self.R
        ep_r = self.ep_r
        n_ns = self.n_ns
        # Compute the fold-change
        pact = MWC(c, ka, ki, ep_ai, n_sites).pact()
        fc = (1 + pact * (R / n_ns) * np.exp(-ep_r)) ** -1
        leakiness = self.leakiness()
        expanded_ka = 1 + c / ka
        expanded_ki = 1 + c / ki

        # Break it into pieces.
        prefactor = -(fc ** 2) * (R / n_ns) * np.exp(-ep_r) * 2 * c * np.exp(-ep_ai)
        numer = (1 / ka) * expanded_ka * expanded_ki ** 2 - (
            1 / ki
        ) * expanded_ka ** 2 * expanded_ki
        denom = (expanded_ka ** 2 + np.exp(-ep_ai) * expanded_ki ** 2) ** 2
        return (2 / (fc - leakiness)) * prefactor * numer / denom

    def compute_properties(self):
        """
        Computes the leakiness, saturation, dynamic range, EC50, and effective hill 
        coefficient for the architecture. Properties are returned as a dictionary. 
        """
        if self.allo == False:
            raise RuntimeError("Available for allosteric molecules only.")

        # Compute the properties.
        leak = self.leakiness()
        sat = self.saturation()
        dyn_rng = self.dynamic_range()
        EC50 = self.ec50()
        Hill = self.effective_hill()

        return {
            "leakiness": leak,
            "saturation": sat,
            "dynamic_range": dyn_rng,
            "EC50": EC50,
            "effective_hill": Hill,
        }

    def bohr_parameter(self):
        r"""
        Computes the Bohr parameter of the form

        bohr = k_BT(log(pact) + log(R / N_ns) + ep_r / k_BT)
        """
        # Compute pact
        if self.allo is True:
            pact = self.mwc.pact()
        else:
            pact = 1
        # Compute and return the Bohr.
        bohr = self.ep_r - np.log(pact) - np.log(self.R / self.n_ns)
        return bohr


def load_constants():
    """Returns a dictionary of various constants incuding binding energies and copy numbers"""
    return dict(
        O1=-15.3,
        O2=-13.9,
        O3=-9.7,
        RBS1147=60,
        RBS446=124,
        RBS1027=260,
        RBS1=1220,
        Nns=4.6e6,
        n_sites=2,
        ep_AI=4.5,
        Ka=139,
        Ki=0.53,
    )
