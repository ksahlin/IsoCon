
# import scipy
# import math
# import random

# def p_correct(k, mu,variance):
#     return 1 - upper_tail_probability_using_CLT(k, mu,variance)

# def upper_tail_probability_using_CLT(k, mu,variance):
#     """
#         Using Lindebergs condition: https://en.wikipedia.org/wiki/Lindeberg%27s_condition
#         For weighted Bernoillis with different p_i

#         upper_tail_probability means the probability of observing a weighted count more or equal to the observed count
#         given that it''s random errors.
#         If this value is really low it means that it''s low probability that
#         base pair is an error --> p(correct) = (1-upper_tail_probability_using_CLT)
#     """
#     z_score = (k - mu)/ math.sqrt(variance)
#     p_value = scipy.stats.norm.sf(z_score) #one-sided
#     return p_value

# def generate_empirical_distr(alpha_j_edit_distance, epsilon, nr_sim_values):
#     empirical_distr = []
#     for i in range(nr_sim_values):
#         simulated_value = 0
#         for x_i in alpha_j_edit_distance:
#             v_i = random.random()
#             if 1.0 - v_i < epsilon[x_i]:
#                 simulated_value += alpha_j_edit_distance[x_i]
#         empirical_distr.append(simulated_value)

#     return sorted(empirical_distr)

# def p_value_quantile_empirical_distribution(alpha_j_edit_distance, epsilon, p_value_threshold, nr_sim_values):
#     p_no_error = reduce(lambda x, y: x*y, map(lambda z: (1- epsilon[z]),epsilon)) 
#     # print("P NO ERROR:", p_no_error)
#     # print(sorted(map(lambda x: alpha_j_edit_distance[x], alpha_j_edit_distance)))
#     # print(sorted(map(lambda x: epsilon[x], epsilon)))

#     approximation_points = []
#     for x_i in alpha_j_edit_distance:
#         approximation_points.append((x_i, alpha_j_edit_distance[x_i]*epsilon[x_i]))
#     sorted_approximation_points =  sorted(approximation_points, key= lambda tup: tup[1], reverse=True)

#     weight_threshold = 0
#     current_p = 1.0
#     current_w = 0

#     plus_mean_weighted_error_of_rest = 0
#     for x_i, w_times_e in sorted_approximation_points:
#         current_p *= epsilon[x_i]
#         if current_p < p_value_threshold:
#             plus_mean_weighted_error_of_rest += w_times_e
#             break
#         else:
#             current_w += alpha_j_edit_distance[x_i]

#         # print("p-val thresh", p_value_threshold, "current_p", current_p, "current_w", current_w)
#     quantile = current_w
#     quantile_new = current_w + plus_mean_weighted_error_of_rest

#     empirical_distr = generate_empirical_distr(alpha_j_edit_distance, epsilon, 5*nr_sim_values)
#     quantile_emp = empirical_distr[-5]
#     # print("quantile:", quantile, "quantile_emp:", quantile_emp)
#     # import matplotlib.pyplot as plt
#     # plt.hist(empirical_distr, 50)
#     # plt.show()
#     return quantile, quantile_new, quantile_emp


# def poisson_cumulative(s, mean):
#     # mean = n*p

#     prob_sum = 0
#     for i in xrange(0, s+1):
#         prob_sum += poisson_probability(i, mean)
#         # prob_sum += mean**i/factorial(i)

#     # return prob_sum*math.e**(-mean)
#     return prob_sum

# def poisson_probability(actual, mean):
#     # naive:   math.exp(-mean) * mean**actual / factorial(actual)

#     # iterative, to keep the components from getting too large or small:
#     p = math.exp(-mean)
#     for i in xrange(actual):
#         p *= mean
#         p /= i+1
#     return p

# def calc_binomial(n_values, p_e):
#     binomial_hashtable = {}  # key is (n,k) tuple and value is probability

#     for n in n_values:
#         # infer what is the largest k we need to
#         # calculate probabilities for set threshold to 0.999
#         # it won't affect the estimations of its 0.999 or 0.9999 for example.
#         k = 0
#         while True:
#             # if n*p_e >= 10:
#             #     # fast normal approximation is good
#             #     pass

#             if p_e <= 1.0/n and n >= 20:
#                 # poisson approximation is good
#                 mu = n*p_e
#                 prob = scipy.stats.poisson.cdf(k, mu)
#                 binomial_hashtable[(n, k)] = prob
#             else:
#                 # exact binomial
#                 prob = scipy.stats.binom.cdf(k, n, p_e)
#                 binomial_hashtable[(n, k)] = prob

#             if prob > 0.999:
#                 break
#             k += 1

#     return binomial_hashtable
