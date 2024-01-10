from tkinter import font
import numpy as np
import matplotlib.pyplot as plt
from sympy import Si

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# 1. Structure of FJC polymer (t=0)
# Parameters
b = 3.0     # Bond length
N = 100     # Number of bonds


# Coordinates arrays
x = np.zeros(N+1)
y = np.zeros(N+1)
z = np.zeros(N+1)

# Generate the random walk
bx = np.random.uniform(-1, 1, 10*N)
by = np.random.uniform(-1, 1, 10*N)
bz = np.random.uniform(-1, 1, 10*N)
normb = np.sqrt(bx**2 + by**2 + bz**2)    # Trajectory norm

# Ignore and normalize the points further than normb
idb = np.where(normb <= 1)[0][0:N]
bx = bx[idb]/normb[idb]
by = by[idb]/normb[idb]
bz = bz[idb]/normb[idb]

# Compute the trajectory
bonds_ini = b * np.array([bx, by, bz])  # Bond vectors
Qini = np.sum(bonds_ini, axis=1)        # Initial Q vector

# 2. Physical Inputs
# kb = 1.38064852e-23   # Boltzmann constant [J/K]
# T = _                 # Temperature [K]
kb_T = 1
# Force vector along x (Fx, 0, 0) [N]
F_array = np.array([[0.01, 0, 0], [0.1, 0, 0], [0.2, 0, 0],
                    [0.35, 0, 0], [0.5, 0, 0], [0.65, 0, 0],
                    [0.82, 0, 0], [1, 0, 0], [1.5, 0, 0],
                    [2, 0, 0], [2.5, 0, 0], [3, 0, 0],
                    [3.5, 0, 0], [4, 0, 0], [4.5, 0, 0],
                    [5, 0, 0], [5.5, 0, 0], [6, 0, 0],
                    [6.5, 0, 0], [7, 0, 0], [7.5, 0, 0],
                    [8, 0, 0], [8.5, 0, 0], [9, 0, 0],
                    [9.5, 0, 0], [10, 0, 0]])

# 8. Extension vs Force analysis
Qx_force = list()      # Final Q depending on force vector
Q_theory = list()     # Theoretical Q depending on force vector
ct = b / kb_T
for ff in F_array:
    alpha = ff[0] * ct  # Force value times constant for alpha
    Q_theory.append(N * b * (1 / np.tanh(alpha) - 1 / alpha))


for F in F_array:
    # 3. Compute initial V
    bonds = bonds_ini.copy()    # Initial bond vectors
    Q0 = Qini.copy()            # Initial Q vector
    V0 = -np.dot(F, Qini)       # Initial pot energy

    # 7. MONTE CARLO PROCESS
    # Iteration parameters
    t_accept = 1    # Time step definition
    naccept = 1000  # Number of accepted steps
    count = 1       # Counter
    nmax = 10000000     # Maximum number of steps
    avg_index = 100     # Average index

    # Sanity check
    if naccept > nmax:
        raise ValueError('naccept, {}, must be less than nmax,\
                          {}'.format(naccept, nmax))

    Q_ext = np.zeros((naccept, 3))  # Q extension vector
    Q_ext[0] = Q0                   # Initial Q extension

    exp_factor = -1 / kb_T   # Pre-calculate this constant

    while t_accept < naccept:
        # 4. Modify bond vector
        bm = np.random.randint(0, N)    # randint is the bond to modify

        # New bond vector
        bx = np.random.uniform(-1, 1, 100)
        by = np.random.uniform(-1, 1, 100)
        bz = np.random.uniform(-1, 1, 100)
        normb = np.sqrt(bx**2 + by**2 + bz**2)    # Trajectory norm
        # Ignore and normalize the points further than normb
        idb = np.where(normb <= 1)[0][0]
        bx = bx[idb]/normb[idb]
        by = by[idb]/normb[idb]
        bz = bz[idb]/normb[idb]

        new_b = b * np.array([bx, by, bz])   # New bond vector

        old_b = bonds[:, bm]    # Old bond vector

        # Compute new Q1 using the new and old bond vectors
        Q1 = Q0 + new_b - old_b

        # 5. Compute new V
        V1 = -np.dot(F, Q1)

        # 6. Acceptance criteria
        delta_V = V1 - V0
        if delta_V <= 0:
            Q0 = Q1
            V0 = V1
            bonds[:, bm] = new_b
            Q_ext[t_accept] = Q1
            t_accept += 1
        # Condition with random number from 0 to 1
        elif np.random.rand() < np.exp(delta_V * exp_factor):
            Q0 = Q1
            V0 = V1
            bonds[:, bm] = new_b
            Q_ext[t_accept] = Q1
            t_accept += 1

        count += 1
        # Sanity check
        if count >= nmax:
            raise ValueError('Maximum number of steps reached, {} accepted \
                              out of {}'.format(t_accept, nmax))

    # Print results
    print('Force: {}'.format(F))
    print('Percentage of acceptance: {} %\n'.format(100 * t_accept/count))
    if F[0] == 10:
        # Plot Q extension
        plt.figure(figsize=(7, 7))
        plt.plot(Q_ext[:, 0], label='Qx')
        plt.title(r'$F = %.2f$' % F[0], fontsize=24)
        plt.xlabel(r'$t$', fontsize=28)
        plt.ylabel(r'$Q_x$', fontsize=28)
        # horizontal line indicating maximum value
        plt.axhline(y=b*N, color='r', linestyle='--',
                    label=r'$Q_x$ Theoretical')
        plt.ylim(0, Q_ext[-1, 0] * 1.2)
        plt.tick_params(axis='y', which='major', labelsize=24, direction='in')
        plt.tick_params(axis='x', which='major', labelsize=24, direction='in')
        plt.show()

    Qx_force.append(np.mean(Q_ext[-avg_index:, 0]))   # Save final Qx value

exten_force = np.array(Qx_force) / (N * b)      # Extension array
exten_theory = np.array(Q_theory) / (N * b)     # Theoretical extension array

# 9. Plot Q vs F
plt.figure(figsize=(7, 7))
plt.plot(exten_force, F_array[:, 0], 'bx', label=r'$Q_x/(n\cdot b)$ Simulated')
plt.plot(exten_theory, F_array[:, 0], 'r--',
         label=r'$Q_x/(n \cdot b)$ Theoretical')
plt.xlabel(r'$Q_x/(n \cdot b)$', fontsize=24)
plt.ylabel(r'$F$', fontsize=24)
plt.tick_params(axis='y', which='major', labelsize=24, direction='in')
plt.tick_params(axis='x', which='major', labelsize=24, direction='in')
plt.legend(fontsize=24, loc='best')
plt.show()
