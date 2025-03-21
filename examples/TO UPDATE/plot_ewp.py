"""
This code plots optical properties.
"""

plt.ion()

#new_dict
plt.subplot(4,1,1)
plt.scatter(op_dict["wavelengths"],op_dict["beta_ext"],label="beta_ext",marker="x")
plt.xlabel("wavelength (um)")
plt.ylabel("beta_ext")

plt.subplot(4,1,2)
plt.scatter(op_dict["wavelengths"],op_dict["ssalb"],label="ssalb",marker="x")
plt.xlabel("wavelength (um)")
plt.ylabel("ssalb")

plt.subplot(4,1,3)
[plt.scatter(np.linspace(0,180,181,endpoint=True),op_dict["phase_function"][i,:],label="phase function",marker="x") for i in [0,int(op_dict["phase_function"].shape[0]/2),op_dict["phase_function"].shape[0]-1]]
plt.xlabel("wavelength (um)")
plt.ylabel("phase_function")
plt.gca().invert_xaxis()


plt.subplot(4,1,4)
[plt.scatter(range(op_dict["legendre_coefficient"].shape[1]),op_dict["legendre_coefficient"][i,:],label="legendre coefficients",marker="x") for i in [0,int(op_dict["legendre_coefficient"].shape[0]/2),op_dict["legendre_coefficient"].shape[0]-1]]
plt.xlabel("wavelength (um)")
plt.ylabel("legendre coefficients")

#old_dict
plt.subplot(4,1,1)
plt.plot(old_dict["wavelengths"],old_dict["beta_ext"],label="beta_ext",c="tab:orange")
plt.xlabel("wavelength (um)")
plt.ylabel("beta_ext")

plt.subplot(4,1,2)
plt.plot(old_dict["wavelengths"],old_dict["ssalb"],label="ssalb",c="tab:orange")
plt.xlabel("wavelength (um)")
plt.ylabel("ssalb")

plt.subplot(4,1,3)
[plt.plot(np.linspace(0,180,181,endpoint=True),old_dict["phase_function"][i,:],label="phase function") for i in [0,int(old_dict["phase_function"].shape[0]/2),old_dict["phase_function"].shape[0]-1]]
plt.xlabel("wavelength (um)")
plt.ylabel("phase_function")


plt.subplot(4,1,4)
[plt.plot(range(old_dict["legendre_coefficient"].shape[1]),old_dict["legendre_coefficient"][i,:],label="legendre coefficients") for i in [0,int(old_dict["legendre_coefficient"].shape[0]/2),old_dict["legendre_coefficient"].shape[0]-1]]
plt.xlabel("wavelength (um)")
plt.ylabel("legendre coefficients")

plt.tight_layout()
plt.show()
