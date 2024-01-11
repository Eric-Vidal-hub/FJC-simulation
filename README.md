# FJC-simulation
The code is developed using _*iPythonNotebook*_ and _*Python*_.
In the first instance, we have to generate multiple configurations which is done in the **file 0 and 1** of an FJC polymer in order to compare with the theoretical results and evaluate the validity of the approximations.
For this purpose, we first study the polymer structure generated, the end-to-end distance, and the gyration radius to verify the model and the simulation are correct.
Once the simulation is done, the mentioned metrics are computed as well as how the error decreases with the number of conformations for $\langle Q^2 \rangle$ and $\langle R_g^2 \rangle$, in **file 2 and 3**.
Moreover, the probability distribution function, $P(Q)$ for $N=100$ bonds, **file 4**, and for the singular behavior $N=2$ bonds, **file 5**, are computed as well as the structure factor $I(k)$ and its Guinier approximiation, in **file 6**
Finally, the Metropolis Monte Carlo algorithm, in **file 7**, is used to simulate the expansion of the polymer when a force is applied along the x-axis, $F_x$.
In this line, the extension in $x$, $Q_x$ is studied in function of the number of steps, considered as time, $t$, and in function of $F_x$.
