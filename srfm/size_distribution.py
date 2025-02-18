"""
Name: size_disitribution.py
Parent package: srfm
Author: Don Grainger
Contributors: Antonin Knizek
Date: 3 January 2025
Purpose: Creates or returns properties of a particle size distribution.
Units: Particle size is expressed in um (micrometers)
       The distribution (n) is in number  per um per cm^3
Integrated values are:
Number density is in number per cm^3
Surface ares density is in um^2 per cm^3
Volume density is in um^3 per cm^3
"""
from abc import ABC, abstractmethod
import numpy as np


# Abstract Base Class
class SizeDistribution(ABC):
    """
    Base class for size distributions.
    """

    def __init__(self, type):
        self.type = type

    @abstractmethod
    def mean(self):
        """
        Calculate the mean of the distribution.
        """
        pass

    # @abstractmethod
    # def other_params(self):
    # """
    # Example method to retrieve other distribution-specific parameters.
    # """
    # pass


# Gaussian Distribution Subclass
class GaussianDistribution(SizeDistribution):
    """
    Creates a Gaussian distribution. NOT IMPLEMENTED
    """

    def __init__(
        self, n=None, r=None, s=None, surface_area_density=None, volume_density=None
    ):
        super().__init__("gaussian")
        self.n = n  # number per cm^3
        self.r = r  # um
        self.s = s

    def mean(self):
        return self.r


# Log-Normal Distribution Subclass
class LogNormalDistribution(SizeDistribution):
    """
    Creates a log-normal distribution.
    Note although the distribution is usually defined by n, r & s the code allows the distribution to be set via
    the surface area density and volume density.
    n - concentration total particles per cm^3
    r - median radius
    s - geometric standard deviation (spread)
    """

    def __init__(
        self, n=None, r=None, s=None, surface_area_density=None, volume_density=None
    ):
        super().__init__("log_normal")

        # Check which combination of inputs is provided
        if n is not None and r is not None and s is not None:
            if n <= 0:
                raise ValueError("n must be greater than 0.")
            if r <= 0:
                raise ValueError("r must be positive.")
            if s <= 0:
                raise ValueError("s must be positive.")
            self.n = n
            self.surface_area_density = (
                n * 4 * np.pi * r**2 * np.exp(2 * np.log(s) ** 2)
            )
            self.volume_density = (
                n * 4 * np.pi * r**3 * np.exp(9 * np.log(s) ** 2 / 2) / 3
            )
        elif surface_area_density is not None and r is not None and s is not None:
            if surface_area_density <= 0:
                raise ValueError("surface_area_density must be greater than 0.")
            if r <= 0:
                raise ValueError("r must be positive.")
            if s <= 0:
                raise ValueError("s must be positive.")
            self.n = surface_area_density / (
                4 * np.pi * r**2 * np.exp(2 * np.log(s) ** 2)
            )
            self.surface_area_density = surface_area_density
            self.volume_density = (
                self.n * 4 * np.pi * r**3 * np.exp(9 * np.log(s) ** 2 / 2) / 3
            )
        elif volume_density is not None and r is not None and s is not None:
            if volume_density <= 0:
                raise ValueError("volume_density must be greater than 0.")
            if r <= 0:
                raise ValueError("r must be positive.")
            if s <= 0:
                raise ValueError("s must be positive.")
            self.n = (
                3 * volume_density / (4 * np.pi * r**3 * np.exp(9 * np.log(s) ** 2 / 2))
            )
            self.surface_area_density = (
                self.n * 4 * np.pi * r**2 * np.exp(2 * np.log(s) ** 2)
            )
            self.volume_density = volume_density
        else:
            # If no valid combination is found, raise an error
            raise ValueError(
                "Invalid parameter combination. Provide one of the following sets:\n"
                "1. n, r, and s\n"
                "2. surface_area_density, r, and s\n"
                "3. volume_density, r, and s."
            )
        self.r = r
        self.s = s
        #       pre-calulate constants to speed up code later
        self.lnr = np.log(r)
        self.lns = np.log(s)
        self.oonls = 1 / self.lns
        self.mootnls = -1 / (2 * self.lns * self.lns)
        self.oostp = 1 / np.sqrt(2 * np.pi)

    def mean(self):
        return self.r * np.exp(0.5 * self.lns ** 2)

    def value(self, radii):
        return (
            self.n
            * self.oostp
            * self.oonls
            * np.exp(self.mootnls * (np.log(radii) - self.lnr) ** 2)
            / radii
        )


#  Selector
def create_distribution(dist_type, **kwargs):
    if dist_type == "log_normal":
        return LogNormalDistribution(**kwargs)
    elif dist_type == "gaussian":
        return GaussianDistribution(**kwargs)
    else:
        raise ValueError(f"Unknown distribution type: {dist_type}")

if __name__ == "__main__":
    # log-normal distribution unit test
    sd = create_distribution("log_normal",
                             n = 1,
                             r = np.exp(1),
                             s = np.exp(1)
                         )
    assert float(sd.value(1)) == 0.24197072451914337, """Log_normal distribution returns
        incorrect value. For n = 1, r = e, s = e, eval. at 1 should be 
        0.24197072451914337."""
    
    sd = create_distribution("log_normal",
                             n = 1,
                             r = np.exp(2),
                             s = np.exp(3)
                         )
    assert float(sd.value(0.5)) == 0.17775476480655455, """Log_normal distribution returns
        incorrect value. For n = 1, r = exp(2), s = exp(3), eval. at 0.5 should be 
        0.17775476480655455."""
