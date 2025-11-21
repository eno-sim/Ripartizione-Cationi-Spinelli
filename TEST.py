import numpy as np
from iminuit import Minuit
from ripcat import RipartoSpinelli
 # Replace 'your_module' with the actual module name where RipartoSpinelli is defined


riparto = RipartoSpinelli.with_defaults()

X0 = np.array([
    0.1253, 0.2006, 0.0264, 0.6451, 0.0020, 0.0000, 0.0000,
    1.3780, 0.0084, 0.0345, 0.1362, 0.0000, 0.4380, 0.0022
])


err0 = np.array([
    0.0020, 0.0049, 0.0039, 0.0036, 0.0010, 0.0000, 0.0000,
    0.0067, 0.0010, 0.0045, 0.0017, 0.0000, 0.0060, 0.0010
])


input_file = "C:\\Users\\enosim\\Desktop\\THESIS\\RIPARTIZIONE-CATIONI-SPINELLI\\SP-3G13.IN"
dis_file = "C:\\Users\\enosim\\Desktop\\THESIS\\RIPARTIZIONE-CATIONI-SPINELLI\\RIP99A.DIS" 

RipartoSpinelli.prepare_input_file(input_file)



j = 0
def fcn_wrapper(*X):
    global j
    j += 1
    X = np.array(X)
    # For the first call, use IFLAG=1 to initialize reading input/dis files
    # (like Fortran CALL FCN(1))
    if not hasattr(riparto, 'initialized'):
        riparto.initialized = True
        F_initial = riparto.FCN(X, 1, input_file, dis_file)
        print(f"Initial F = {F_initial}")
        print('iteration ', j)
        return F_initial
    else:
        # Subsequent calls use IFLAG=3 (as in MIGRAD / IMPROVE / CALL FCN 3)
        print('iteration ', j)
        return riparto.FCN(X, 3, input_file, dis_file)


m = Minuit(fcn_wrapper, *X0)

# Set limits or constraints if needed
for i in range(len(X0)):
    if i < len(X0)//2:
        m.limits[i] = (0,1)  # example bounds; tune as needed
    else:
        m.limits[i] = (0,2)


# Optionally set strategy, tolerance, etc.
m.strategy = 2
m.tol = 1e-4
m.errordef = 1e-4
m.simplex()
m.migrad()



# CALL FCN 3 equivalent — evaluate FCN with the best-fit parameters
final_F = riparto.FCN(np.array(m.values), 3, input_file, dis_file)
print(f"Final minimized F = {final_F}")


print("Best-fit parameters:")
print(m.values)
print("Parameter errors:")
print(m.errors)
print('Limits:')
print(m.limits)
