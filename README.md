# Computational-Systems-Biology-2020
Semester course project- Computational Systems Biology (BT5240)

The project is a study of an automaton model of the cell cycle coupled with a basic biochemical ODE model that progresses through the phases: G1, S, G2 and M. The transition between each phase controlled by the concentrations of the major molecules, total cyclin and MPF, which in turn are regulated by other molecules.
Cells on completing M phase, divide into two daughter cells which enter into the G1 phase of the cycle. Further, during the cycle each cell has a defined probability for exiting the cycle and undergo apoptosis. The exit can take place only during G1-S and G2-M transitions thus establishing the role of molecular check points in a cycle.
To achieve homeostasis and the maintenance of total number of cells, we assumed a fixed cell population with apoptosis balancing cell replication. We also assumed that the cell divides symmetrically and has not considered cell size variations.
With this crude model, we attempted to study the dynamics of the cell cycle such as the effect of the concentrations of the molecules, the number and the evolution of cell cycle phases in the cell population with time, mean duration of phases, number of cells in the population, and their steady-state behaviours corresponding to the different initial conditions.
