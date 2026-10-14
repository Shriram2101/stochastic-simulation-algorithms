
# Exact Stochastic Simulation of Coupled Chemical Reactions

This project implements the stochastic simulation algorithm introduced by Daniel T. Gillespie in his paper *Exact Stochastic Simulation of Coupled Chemical Reactions* (1977). The objective is to reproduce selected numerical results presented in the paper and explore the application of exact stochastic simulation to chemical reaction systems.

## 1. Overview of the Paper

### Abstract Summary

Chemical reaction systems are often modeled using deterministic reaction-rate equations, which describe the evolution of species concentrations through coupled ordinary differential equations. Although this approach is useful for many macroscopic systems, it does not explicitly account for the discrete and random nature of molecular interactions.

In his 1977 paper, Daniel T. Gillespie introduced a stochastic simulation algorithm that numerically generates the time evolution of a well-mixed chemical system using a Monte Carlo procedure. Rather than solving the chemical master equation directly, the algorithm determines the time of the next reaction event and selects which reaction occurs based on the reaction propensities.

The method captures molecular population fluctuations and correlations that deterministic descriptions may overlook. It also advances time by sampling the intervals between reaction events rather than approximating continuous time with fixed, finite time steps. Gillespie demonstrated the method using established model systems, including the Lotka model, the Brusselator, and the Oregonator.

This project aims to implement the algorithm and reproduce selected results from the original paper.

## 2. Gillespie's Stochastic Simulation Algorithm

The Gillespie algorithm, also known as the stochastic simulation algorithm (SSA), is an event-driven method for simulating the time evolution of a well-mixed system of chemical reactions.

Consider a system containing $N$ chemical species participating in $M$ reaction channels. The state of the system at time $t$ is represented by the molecular population vector

$$
\mathbf{X}(t) = (X_1(t), X_2(t), \ldots, X_N(t)),
$$

where $X_i(t)$ is the number of molecules of species $i$.

Each reaction channel $\mu$ is associated with a propensity function $a_\mu(\mathbf{X})$. Given the current state, $a_\mu(\mathbf{X})\,dt$ gives the probability, to first order in $dt$, that reaction $\mu$ occurs during the next infinitesimal time interval $dt$.

### Algorithm

At each iteration, the algorithm performs the following steps:

1. **Initialize the system.** Specify the initial molecular populations, reaction channels, and kinetic parameters. Set the initial time to $t=0$.

2. **Calculate the reaction propensities.** Evaluate the propensity $a_\mu(\mathbf{X})$ for each reaction channel and calculate the total propensity
   $a_0(\mathbf{X}) = \sum_{\mu=1}^{M} a_\mu(\mathbf{X})$

3. **Generate random numbers.** Draw two independent uniform random numbers $r_1,r_2\in(0,1)$.

4. **Determine the time to the next reaction.** Sample the waiting time from an exponential distribution:\
   $\tau = \frac{1}{a_0(\mathbf{X})}\ln\left(\frac{1}{r_1}\right)$

5. **Select the reaction channel.** Choose reaction $\mu$ such that\
   $\sum_{\nu=1}^{\mu-1} a_\nu(\mathbf{X})< r_2a_0(\mathbf{X})\leq\sum_{\nu=1}^{\mu} a_\nu(\mathbf{X})$

7. **Update the molecular populations.** Execute the selected reaction, update the state vector according to the reaction's stoichiometry, and advance time by $t\leftarrow t+\tau$.

8. **Repeat.** Recalculate the propensities and continue until the specified simulation stopping condition is reached.

The algorithm generates a statistically exact realization of the specified well-mixed, continuous-time Markov reaction model, subject to the validity of that model and the accuracy of its implementation.

### Why use the Gillespie algorithm?

- **Discrete molecular populations:** Molecular counts change through individual reaction events.
- **Intrinsic stochasticity:** Randomness in reaction events is explicitly represented.
- **Exact event timing:** Reaction waiting times are sampled from the model's probability distribution rather than approximated using a fixed time step.
- **Fluctuations and correlations:** Individual trajectories reveal stochastic behavior that may not be apparent from deterministic reaction-rate equations.

## 3. Reproducing the Results

The main objective of this project is to reproduce selected numerical results from Gillespie's original paper.

The implementation and analysis will include:

- Defining the reaction networks and kinetic parameters.
- Implementing the stochastic simulation algorithm.
- Simulating the relevant model systems.
- Comparing simulated trajectories and statistical behavior with the published results.

## 4. Reference

Gillespie, D. T. (1977). Exact Stochastic Simulation of Coupled Chemical Reactions. *The Journal of Physical Chemistry*, **81**, 2340.

The original publication is the primary reference for the algorithm and the numerical results investigated in this project.
