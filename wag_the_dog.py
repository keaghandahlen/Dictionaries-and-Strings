import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def wag_the_dog(potential_name, adjustable_coefficient_values, length_scale, positions, initial_values):
    # List of currently supported potentials
    supported_potentials = ['infinite square well']

    # Set up tuple for the range of positions to use in solve_ivp
    position_range = (np.min(positions), np.max(positions))

    if potential_name == 'infinite square well':
        from infinite_square_well import infinite_square_well_differential_equation, draw_infinite_square_well_potential
        for adjustable_coefficient in adjustable_coefficient_values:
            # Assign K to an list
            arguments = [adjustable_coefficient]
            # Solve y''(x) = f[x, y(x), y'(x), K] given (x_min, x_max), (y(x_min), y'(x_min)), x, and K
            solution = solve_ivp(infinite_square_well_differential_equation, position_range, initial_values,
                                 t_eval=positions, args=arguments)
            # Assign the normalized solution to psi to the y-axis
            #   A = 1 / sqrt( int( y(x)^2 ) dx )
            normalization_factor = 1. / np.sqrt(np.trapz(solution.y[0]**2, x=positions))
            #   psi(x) = A * y(x)
            psi_numerical = normalization_factor * solution.y[0]
            # Plot normalized solution
            plt.plot(positions, psi_numerical, label=r'$k^2 = {:.2f}$'.format(adjustable_coefficient))

        plt.legend()
        draw_infinite_square_well_potential(positions, infinite_square_well_width=length_scale)

    else:
        print('No potential named {} supported'.format(potential_name))
        print('Try one from this list: {}'.format([potential for potential in supported_potentials]))
        return

    # Formatting for all plots
    #    Draw psi = 0 line
    plt.axhline(color='black')
    #   Label axes
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\psi_{n}(x)$')
    #   Show plot
    plt.show()

    return
