def calculate_age_distribution_properties():
    import matplotlib.pyplot as plt
    import numpy as np

    age_distribution = {
        14: 1,
        15: 1,
        16: 3,
        22: 2,
        24: 2,
        25: 5,
    }

    # j
    ages = [*age_distribution]
    # N(j)
    numbers = [age_distribution[age] for age in ages]

    # Plot N(j) vs j
    plt.bar(ages, numbers)
    plt.ylim(0, 6)
    plt.ylabel(r'$N(j)$')
    plt.xlim(10, 27)
    plt.xlabel(r'$j$')
    plt.show()

    # N
    total_number = 0.
    print('N = {}'.format(total_number))

    # P(15)
    probability_distribution = np.zeros(25)
    print('P(15) = {}'.format(probability_distribution[15]))

    # sum of probabilities
    probability_sum = 0.
    print('sum(P) = {}'.format(probability_sum))

    # most probable age
    most_probable_age = 0.
    print('max(P) = {}'.format(most_probable_age))

    # median age
    age = 0
    print('P reaches 0.5 at j={}'.format(age))

    # average age
    mean_age = 0.
    # or:
    mean_age2 = 0.
    print('<j> = {} = {}'.format(mean_age, mean_age2))

    # expectation value of j^2
    age_squared_expectation = 0.
    print('<j^2> = {}'.format(age_squared_expectation))
    print('<j>^2 = {}'.format(mean_age ** 2))

    # Δj
    delta_j = 0.
    print('Δj = {}'.format(delta_j))

    # <Δj>
    delta_j_expectation = 0.
    print('<Δj> = {}'.format(delta_j_expectation))

    # σ^2
    probability_distribution_array = 0.
    variance = 0.
    print('σ^2 = {} = {}'.format(variance, age_squared_expectation - mean_age ** 2))

    # σ
    standard_deviation = 0.
    print('σ = {}'.format(standard_deviation))

    return
