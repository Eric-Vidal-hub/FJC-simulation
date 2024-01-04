import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# 1. Structure of FJC polymer (t=0)
# Parameters
b = 3.0     # Bond length
N = 100     # Number of bonds
T = 100     # Number of conformations

# Coordinates arrays
x = np.zeros((T, N+1))
y = np.zeros((T, N+1))
z = np.zeros((T, N+1))

# Read simulation data
filename = 'simulation_FJC_b=%.1f_N=%d_T=%d.xyz' % (b, N, T)
with open(filename, 'r', encoding='utf-8') as f:
    for t in range(T):
        # First 2 lines unnecessary
        lines = f.readline()
        lines = f.readline()
        # Save coordinates and separate in x, y, z
        for i in range(N+1):
            lines = f.readline()
            coord = lines.split()
            x[t, i] = float(coord[1])
            y[t, i] = float(coord[2])
            z[t, i] = float(coord[3])

# Initial polymer for one T (t=0)
polymer = np.array([x[0], y[0], z[0]])
# Bonds vector
bonds = np.array(polymer[:, 1:] - polymer[:, :-1])

# 2. Physical Inputs
kb = 1.38064852e-23  # Boltzmann constant [J/K]
T = 1e21              # Temperature [K]
kb_T = kb * T
# Force vector along x (Fx, 0, 0) [N]
# F_array = np.array([[0.01, 0, 0], [0.1, 0, 0], [0.5, 0, 0],
#                     [1, 0, 0], [1.5, 0, 0], [2, 0, 0],
#                     [2.5, 0, 0], [3, 0, 0], [3.5, 0, 0],
#                     [4, 0, 0], [4.5, 0, 0], [5, 0, 0],
#                     [5.5, 0, 0], [6, 0, 0], [6.5, 0, 0],
#                     [7, 0, 0], [7.5, 0, 0], [8, 0, 0],
#                     [8.5, 0, 0], [9, 0, 0], [9.5, 0, 0],
#                     [10, 0, 0]])

# F_array = np.array([[0.01, 0, 0], [0.1, 0, 0], [0.5, 0, 0],
#                     [1, 0, 0], [2, 0, 0], [3, 0, 0]])

F_array = np.array([[0.01, 0, 0], [0.1, 0, 0], [1, 0, 0]])

# # Same left plot as him
# F_array = np.array([[0.1, 0, 0]])

# 8. Extension vs Force analysis
Qx_force = list()      # Final Q depending on force vector
Q_theory = list()     # Theoretical Q depending on force vector
ct = b / kb_T
for ff in F_array:
    alpha = ff[0] * ct  # Force value times constant for alpha
    Q_theory.append(N * b * (1 / np.tanh(alpha) - 1 / alpha))


for F in F_array:
    # 3. Compute initial V
    Q0 = np.array(polymer[:, -1] - polymer[:, 0])
    V0 = -np.dot(F, Q0)

    # 7. MONTE CARLO PROCESS
    # Iteration parameters
    t_accept = 1    # Time step definition
    naccept = 1000  # Number of accepted steps
    count = 1       # Counter
    nmax = 10000000    # Maximum number of steps

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
        new_b = np.random.rand(3)   # Generate 3 random numbers at once
        new_b *= b / np.linalg.norm(new_b)  # Normalize and scale in-place

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
    # if F[0] == 0.01:
    # Plot Q extension
    plt.figure()
    plt.plot(Q_ext[:, 0], label='Qx')
    plt.title(r'$F = %.2f$' % F[0], fontsize=18)
    plt.xlabel(r'$t$', fontsize=18)
    plt.ylabel(r'$Q_x$', fontsize=18)
    plt.tick_params(axis='y', which='major', labelsize=16, direction='in')
    plt.tick_params(axis='x', which='major', labelsize=16, direction='in')
    plt.show()

    Qx_force.append(Q_ext[-1, 0])   # Save final Qx value

# 9. Plot Q vs F
plt.figure()
plt.plot(F_array[:, 0], Qx_force, 'bo', label=r'$Q_x$ Simulated')
plt.plot(F_array[:, 0], Q_theory, 'r--', label=r'$Q_x$ Theoretical')
plt.xlabel(r'$F$', fontsize=18)
plt.ylabel(r'$Q_x$', fontsize=18)
plt.tick_params(axis='y', which='major', labelsize=16, direction='in')
plt.tick_params(axis='x', which='major', labelsize=16, direction='in')
plt.legend()
plt.show()
