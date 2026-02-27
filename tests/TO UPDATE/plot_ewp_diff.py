"""
This code snippet attempts to plot differences in optical properties (plots the op_dict)
"""

plt.figure(figsize=(11.7, 8.3))

plt.ion()
plt.subplot(3, 1, 1)
plt.plot(diff_dict["wavelengths"], diff_dict["beta_ext"], label="beta_ext")
plt.xlabel("wavelength (um)", fontsize="12")
plt.ylabel(r"$\Delta \beta ^{ext} (\lambda)$ (%)", fontsize="12")
# plt.yscale("log")
plt.title("Extinction coefficient", fontsize="15")

plt.subplot(3, 1, 2)
plt.plot(diff_dict["wavelengths"], diff_dict["ssalb"], label="ssalb")
plt.xlabel("wavelength (um)", fontsize="12")
plt.ylabel(r"$\Delta \omega (\lambda)$ (%)", fontsize="12")
# plt.yscale("log")
plt.title("Single scattering albedo", fontsize="15")


plt.subplot(3, 1, 3)
[
    plt.plot(
        diff_dict["wavelengths"],
        diff_dict["phase_function"][:, i],
        label=f"phase function, {i} deg",
    )
    for i in range(diff_dict["phase_function"].shape[1])
]
plt.xlabel("wavelength (um)", fontsize="12")
plt.ylabel(r"$\Delta p(\mu, \lambda)$ (%)", fontsize="12")
# plt.yscale("log")
# plt.legend()
plt.title("Phase function - each line represents one angle", fontsize="15")

# plt.subplot(4,1,4)
# [plt.plot(diff_dict["wavelengths"],diff_dict["legendre_coefficient"][:,i],label=f"legendre coefficient no. {i}") for i in range(diff_dict["legendre_coefficient"].shape[1])]
# plt.xlabel("wavelength (um)")
# plt.ylabel("legendre coefficients")
# plt.yscale("log")
##plt.legend()


plt.suptitle(
    "Difference in optical properties arising from interpolation", fontsize="18"
)
plt.tight_layout()
plt.show()
# plt.legend()
