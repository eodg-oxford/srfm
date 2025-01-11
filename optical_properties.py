# particle_optical_properties.py

import numpy as np
import quadrature as quad
import ARIA_module as ARIA  # Import the RI class from ri_module
import mie_module  # Assuming mie_module is the compiled Fortran module

def legendre_polynomial_expansion(inp, qv, qw, phase):
    """
    Expands a function as a Legendre series. The function is assumed to have been
    evaluated at Legendre quadrature points. The evaluation stops when the number 
    of terms exceeds the number of quadrature points or when the absolute value of 
    the Legendre coefficient is less than 1e-9.

    Parameters:
        inp (int): Number of points.
        qv (array-like): Legendre points.
        qw (array-like): Legendre weights.
        phase (array-like): Input (phase) function.

    Returns:
        tuple: lc (Legendre coefficients), inlc (Number of coefficients used).
    """
    Imaxnp = 20000
    if inp > Imaxnp:
        raise ValueError("Error in legendre_polynomial_expansion: Too many quadrature points")

    # Initialize Legendre coefficients
    lc = np.zeros(inp, dtype=np.float64)

    # Compute the first two Legendre coefficients
    lc[0] = np.sum(phase * qw) / 2.0
    lc[1] = 1.5 * np.sum(phase * qv * qw) / 2.0

    # Initialize Legendre polynomials
    lpnm2 = np.ones(inp, dtype=np.float64)  # P_0 = 1
    lpnm1 = qv                             # P_1 = qv

    # Iterate to compute higher-order coefficients
    n = 2
    while n < inp:
        # Compute the nth Legendre polynomial
        lpn = ((2 * n - 1) / n) * qv * lpnm1 - ((n - 1) / n) * lpnm2
        lpnm2 = lpnm1
        lpnm1 = lpn

        # Compute the nth Legendre coefficient
        lc[n] = (2 * n + 1) * np.sum(phase * lpn * qw) / 2.0

        # Stop if the coefficient is below the threshold
        if abs(lc[n]) < 1e-9:
            break

        n += 1

    # Number of coefficients used
    inlc = n - 1
    return lc, inlc



# Determine the optical properties of a distribution of particles as a function of composition and structure
def ewp_hs(wavelength, composition, distribution, refractive_index=None, angle = None, legendre_coefficients_flag=False, radii =200, eta = 1e-6):
    # return the extinction, singlescatter albedo and phase function for a particle size of a given 'distribution' and 'composition'
    # the phase function is given at the specified angles if given else at 1 degree resolution
    # if the legendre_coefficients_flag is True angle is ignored and Lobatto quadrature angles are used

    wavelengths = wavelength.size
  
    
    # If angles is given calculate the phase function at those angles otherwise use default 0 to 180 in steps of 1 degree
    if legendre_coefficients_flag:
       cos_angle_value, cos_angle_weight  = quad.quadrature101("L", 181)    # change quadrature N if necessary         
       if angle is not None: 
           print('Warning in ewp_hs: Specified angles ignored as Legendre coefficients requested')          
    else:
       if angle is not None: 
           cos_angle_value = np.cos(angle*np.pi/180)
           cos_angle_weight = np.zeros(cos_angle_value.size)
       else:           
           cos_angle_value, cos_angle_weight  = (0,1) # 1 angle
       
    cos_angle_value =  np.array(cos_angle_value, dtype=np.float64) # ensure correct type for Mie call
    angles = cos_angle_value.size    

    radius_lower_bound = np.exp(np.log(distribution.r) - np.log(distribution.s)**2 - np.sqrt(-2 * np.log(eta) * np.log(distribution.s)**2))
    radius_upper_bound = np.exp(np.log(distribution.r) - np.log(distribution.s)**2 + np.sqrt(-2 * np.log(eta) * np.log(distribution.s)**2))


    radius, weight_temp = quad.quadrature("T", radii, radius_lower_bound, radius_upper_bound)            
    radius_weight = weight_temp * distribution.value(radius)

    if composition == "ri": 
        if refractive_index is None:
            print('Error: refractive index not defned')
        else:
            ri_n = np.array(np.full(wavelengths, np.real(refractive_index)))
            ri_k = np.array(np.full(wavelengths, np.imag(refractive_index)))                   
    else:
        ri_object = ARIA.RI()
        ri_n, ri_k = ri_object.load_refractive_indices(composition, wave=wavelength, mode="wavelength")
        ri_k=-ri_k  
    # Validate refractive index before proceeding
    if np.any(ri_k > 0):
        raise ValueError("Refractive index (k) contains invalid values. Ensure all values are zero or negative.")

    # Initialize arrays with the same size as 'wavelengths'
    beta_ext = np.zeros(wavelengths)  # 1D array for extinction
    beta_sca = np.zeros(wavelengths)  # 1D array for scattering
    phase_function = np.zeros((wavelengths, angles))  # 2D array for phase function

    # Initialize arrays with the same size as 'radii'
    Q_ext_value = np.zeros(radii)  # 1D array for extinction per radius
    Q_sca_value = np.zeros(radii)  # 1D array for scattering per radius
    phase_function_value = np.zeros((radii, angles))  # 2D array for phase function per radius   
   
    Q_ext = np.zeros(1, dtype=np.float64)  # Single output scalar wrapped in an array
    Q_sca = np.zeros(1, dtype=np.float64)  # Single output scalar wrapped in an array
    phase_function_particle = np.zeros(angles, dtype=np.float64)  # Complex output array
    Error = np.array([0], dtype=np.int32)  # Error code as an integer array   

    if legendre_coefficients_flag:
        legendre_coefficient = np.zeros((wavelengths,10000))
 
    for i in range(wavelengths):
        refractive_index = complex(ri_n[i], ri_k[i])
            
        for j in range(radii):        
            size_parameter = 2 * np.pi * radius[j] / wavelength[i]

            # Perform Mie calculations
            Q_ext, Q_sca, phase_function_particle, error = mie_module.mie_ewp(size_parameter, refractive_index, cos_angle_value)
            
            # Store results for this radius
            Q_ext_value[j] = Q_ext
            Q_sca_value[j] = Q_sca
            phase_function_value[j, :] = phase_function_particle  # Ensure correct slicing
        
        # Weighted sums for extinction and scattering
        beta_ext[i] = np.sum(radius_weight * np.pi * radius**2 * Q_ext_value) * 1E-6 # 1/m
        beta_sca[i] = np.sum(radius_weight * np.pi * radius**2 * Q_sca_value) * 1E-6 # 1/m
         
        # Weighted sums for phase function
        for k in range(angles):
            phase_function[i, k] = np.sum(radius_weight * phase_function_value[:, k])
              
        # If requested expand phase functions as Legendre coefficients    
        if legendre_coefficients_flag:
           legendre_coefficient_temp, legendre_coefficient_number = legendre_polynomial_expansion(angles, cos_angle_value, cos_angle_weight, phase_function[i,:])
           legendre_coefficient[i,0:legendre_coefficient_number] = legendre_coefficient_temp[0:legendre_coefficient_number]

    # Return the computed results
    if legendre_coefficients_flag:
        return beta_ext, beta_sca / beta_ext, phase_function, legendre_coefficient
    else:
        return beta_ext, beta_sca / beta_ext, phase_function
