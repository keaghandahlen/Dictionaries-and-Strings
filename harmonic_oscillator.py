import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial, hermite


def calculate_harmonic_oscillator_eigenfunctions(eigenfunction_positions, eigenfunction_index,
                                                 harmonic_oscillator_width=1.0):
    """
    :param eigenfunction_positions:  x-values in a Numpy array
    :param eigenfunction_index: n
    :param harmonic_oscillator_width: a
    :return: psi_n(x) from x_minimum to x_maximum
    """
    # Assuming k = 4 / (2p) [2p = latus rectum of parabola], ω^2 = k/m (m = 1)
    force_constant = 1. / harmonic_oscillator_width
    angular_frequency = np.sqrt(force_constant)

    # A_n = (m ω / π ħ)^(1/4) (1 / sqrt(2^n n!)), m = ħ = 1
    normalization_factor = np.power(angular_frequency / np.pi, 1 / 4)
    normalization_factor /= np.sqrt(np.power(2., eigenfunction_index) * factorial(eigenfunction_index))
    print('A = {}'.format(normalization_factor))

    # ξ = sqrt(m ω / ħ) x, m = ħ = 1
    xi = np.sqrt(angular_frequency) * eigenfunction_positions

    # psi_n(x) = A_n H_n(xi) exp(-xi^2 / 2)
    harmonic_oscillator_eigenfunction_values = normalization_factor * hermite(eigenfunction_index)(xi)
    harmonic_oscillator_eigenfunction_values *= np.exp(-0.5 * xi ** 2)

    return harmonic_oscillator_eigenfunction_values


def draw_harmonic_oscillator_potential(function_positions, harmonic_oscillator_width=1.0,
                                       maximum_potential_value=2.0, plot_axes=np.array([])):
    """
    :param plot_axes:
    :param function_positions:  x-values in a Numpy array
    :param harmonic_oscillator_width:  a
    :param maximum_potential_value: sets maximum V (and minimum V = -maximum V)
    :return:
    """
    force_constant = 1. / harmonic_oscillator_width
    minimum_position = np.min(function_positions)
    minimum_potential_position = minimum_position - 0.333 * harmonic_oscillator_width
    maximum_position = np.max(function_positions)
    maximum_potential_position = maximum_position + 0.333 * harmonic_oscillator_width

    potential_positions = np.linspace(minimum_potential_position, maximum_potential_position)
    potential_values = (1 / 2) * force_constant * potential_positions ** 2
    if plot_axes.size != 0:
        number_of_subplots = np.sum(plot_axes.shape)
        for axis_index in np.arange(number_of_subplots):
            plot_axes[axis_index].fill_between(potential_positions, potential_values,
                                               facecolor='gray', alpha=0.2)
            plot_axes[axis_index].fill_between(potential_positions,
                                               np.full(len(potential_positions), -maximum_potential_value),
                                               facecolor='gray', alpha=0.2)
            plot_axes[axis_index].set_xlim([np.min(potential_positions), np.max(potential_positions)])
            plot_axes[axis_index].set_ylim([-maximum_potential_value, maximum_potential_value])

    else:
        plt.fill_between(potential_positions, potential_values,
                         facecolor='gray', alpha=0.2)
        plt.fill_between(potential_positions, np.full(len(potential_positions), -maximum_potential_value),
                         facecolor='gray', alpha=0.2)

        plt.xlim([np.min(potential_positions), np.max(potential_positions)])
        plt.ylim([-maximum_potential_value, maximum_potential_value])

    return


def harmonic_oscillator_differential_equation(xi, psi, k):
    """
    Sets up second-order differential equation for solve_ivp as two linear differential equations
    psi' = psi'
    psi'' = (xi**2 - K)*psi
    :param xi:  position value
    :param psi:  list of wave function and wave-function derivative's values
    :param K: adjustable parameter to find solutions
    :return: psi': list of wave function and wave-function derivative's first derivative values
    """
    dpsi = [0, 0]  # initialize the (psi', psi'') vector
    dpsi[0] = psi[1]  # set psi' = psi'
    dpsi[1] = (xi ** 2 - k) * psi[0]  # set psi'' = (xi**2 - constant) * psi
    return dpsi
