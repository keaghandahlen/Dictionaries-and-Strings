import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler


def plot_wave_function(potential_name, length_scale, positions, quantum_numbers, plot_probability=False):
    """
    Plots one-dimensional wave function given name of potential, principal length scale of potential, and
    quantum numbers for which to plot
    :param potential_name:
    :param length_scale:
    :param positions:
    :param quantum_numbers:
    :param plot_probability:
    :return:
    """

    # Set curves to cycle through the following colors and styles
    curve_colors = ['blue', 'orange', 'green', 'red']
    curve_styles = ['solid', 'dashed', 'dashdot', 'dotted']
    style_cycler = (cycler('color', curve_colors) +
                    cycler('linestyle', curve_styles))

    # List of currently supported potentials
    supported_potentials = ['infinite square well', 'harmonic oscillator']

    # Go through each potential
    if potential_name == supported_potentials[0]:
        from infinite_square_well import calculate_infinite_square_well_eigenfunctions, \
            draw_infinite_square_well_potential
        if plot_probability:
            figure, axes = plt.subplots(1, 2)
            for quantum_number in quantum_numbers:
                function_values = calculate_infinite_square_well_eigenfunctions(positions, quantum_number,
                                                                                infinite_square_well_width=length_scale)
                # Cycle styles
                # plt.rc('axes', prop_cycle=style_cycler)
                axes[0].set_prop_cycle(style_cycler)
                axes[1].set_prop_cycle(style_cycler)
                # Plot wave functions
                axes[0].plot(positions, function_values)
                # Plot probability densities
                axes[1].plot(positions, function_values**2)

        else:
            for quantum_number in quantum_numbers:
                function_values = calculate_infinite_square_well_eigenfunctions(positions, quantum_number,
                                                                                infinite_square_well_width=length_scale)
                plt.rc('axes', prop_cycle=style_cycler)
                plt.plot(positions, function_values)

        draw_infinite_square_well_potential(positions, infinite_square_well_width=length_scale,
                                            maximum_potential_value=2. / length_scale)

    elif potential_name == supported_potentials[1]:
        from harmonic_oscillator import calculate_harmonic_oscillator_eigenfunctions, \
            draw_harmonic_oscillator_potential
        plot_maximum = 0
        if plot_probability:
            figure, axes = plt.subplots(1, 2)
            # Cycle styles
            axes[0].set_prop_cycle(style_cycler)
            axes[1].set_prop_cycle(style_cycler)

            for quantum_number in quantum_numbers:
                function_values = calculate_harmonic_oscillator_eigenfunctions(positions, quantum_number,
                                                                               harmonic_oscillator_width=length_scale)
                # plt.rc('axes', prop_cycle=style_cycler)
                # Plot wave functions
                axes[0].plot(positions, function_values)
                # Plot probability densities
                axes[1].plot(positions, function_values**2)
                wave_function_maximum = np.amax(function_values)
                plot_maximum = max(wave_function_maximum, plot_maximum)

            plot_maximum *= 1.2
            draw_harmonic_oscillator_potential(positions, harmonic_oscillator_width=length_scale,
                                                   maximum_potential_value=plot_maximum, plot_axes=axes)
            axes[0].set_aspect(1)
            axes[1].set_aspect(1)
        else:
            for quantum_number in quantum_numbers:
                function_values = calculate_harmonic_oscillator_eigenfunctions(positions, quantum_number,
                                                                               harmonic_oscillator_width=length_scale)
                wave_function_maximum = np.amax(function_values)
                plot_maximum = max(wave_function_maximum, plot_maximum)

                plt.rc('axes', prop_cycle=style_cycler)
                plt.plot(positions, function_values)

            plot_maximum *= 1.2
            draw_harmonic_oscillator_potential(positions, harmonic_oscillator_width=length_scale,
                                                   maximum_potential_value=plot_maximum)

    else:
        print('No potential named {} supported'.format(potential_name))
        print('Try one from this list: {}'.format([potential for potential in supported_potentials]))
        return

    # Formatting for all potentials and wave functions
    if plot_probability:
        #    Draw psi = 0 line
        axes[0].axhline(color='black')
        axes[1].axhline(color='black')
        #   Label axes
        axes[0].set_xlabel(r'$x$')
        axes[1].set_xlabel(r'$x$')
        axes[0].set_ylabel(r'$\psi_{n}(x)$')
        axes[1].set_ylabel(r'$|\psi_{n}(x)|^2$')
        #   Plot y-axis ticks and labels on right for subplot 2
        axes[1].yaxis.set_label_position("right")
        axes[1].yaxis.tick_right()
    else:
        #    Draw psi = 0 line
        plt.axhline(color='black')
        #   Label axes
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\psi_{n}(x)$')

    figure.tight_layout()
    #   Make legend
    plt.legend([r'$n = $' + str(quantum_number) for quantum_number in quantum_numbers])
    #   Display figure
    plt.show()

    return
