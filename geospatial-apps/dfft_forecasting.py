from numpy.linalg import matrix_power
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
"""
We shall do a one-time optimization over the change in the chemical potential, mu, and the
time scale of matrix multiplications per year, R, such that the forecast for each county of the 2010
data from the year 2000 and the known change in the county composition is optimized.

To do so, we impose the following constraint on T and mu...
1. The change in composition is fixed. For example, if a city has 30% White in 2000 and 25% in 2010, then
we fix that the forecast must also match a decrease of 5%.
2. The number of times the transition matrix is applied should be the same across all neighborhoods

To search, we shall track the log-likelihood and the expected countywide fraction as we apply
the transition matrix more times. For a given mu, we stop when we get to the correct change in fraction
and report back the final log-likelihood. This log-likelihood is what we optimize to find the appropriate
time scale and chemical potential.

In principle, to forecast into the future, we can assume the time scale stays constant and then find the
correct chemical potential that gives the desired change in fraction.

"""


def estimate_mu(race,df_blockgroup_county,tau,delta_years,hbars,n_county_final,mu_lims = [-1, 1],max_iterations = 10,tol = 0.01):
    # TODO: Add in some error checking, return a value if mu does not converge
    total_mean = int(np.round(df_blockgroup_county["total_10"].mean()))
    Ns = df_blockgroup_county[race+"_10"]
    totals = df_blockgroup_county["total_10"]
    state_vector = np.zeros(total_mean+1)
    total_county = sum(totals)
    for (N,total) in zip(Ns,totals):
        N_transformed = int(round(N*total_mean/total))
        state_vector[N_transformed] += total/total_county
    num_multiplications = int(round(total_mean*tau*delta_years))
    ns_discrete = np.linspace(0,1,total_mean+1)
    # print(hbars)
    H_discrete = interpolate_H_discrete_from_hbar(hbars[hbars["race"] == race]["hbar"].iloc[0],
                                                  hbars[hbars["race"] == race]["ns"].iloc[0], total_mean)
    H_discrete_mu = H_discrete + np.mean(mu_lims)*ns_discrete
    transition_matrix = build_transition_matrix(H_discrete_mu, total_mean, option="right_multiply")
    final_state = np.matmul(np.linalg.matrix_power(transition_matrix,num_multiplications),state_vector)
    final_fraction = np.sum(np.matmul(ns_discrete,final_state))
    iter = 1
    # print("iter: ", iter, " final_fraction: ", final_fraction, " mu: ", np.mean(mu_lims), " target: ", n_county_final)
    # append!(progress_df,DataFrame(mu=mean(mu_lims), final_fraction = final_fraction))
    while (abs(n_county_final-final_fraction)>tol) & (iter<max_iterations):
        iter+=1
        if final_fraction<n_county_final:
            mu_lims[1] = np.mean(mu_lims) #Decrease mu
        else:
            mu_lims[0] = np.mean(mu_lims) #Increase mu
        H_discrete_mu = H_discrete + np.mean(mu_lims) * ns_discrete
        transition_matrix = build_transition_matrix(H_discrete_mu, total_mean, option="right_multiply")
        final_state = np.matmul(np.linalg.matrix_power(transition_matrix,num_multiplications),state_vector)
        final_fraction = np.sum(np.matmul(ns_discrete,final_state))
        # print("iter: ", iter, " final_fraction: ", final_fraction, " mu: ", np.mean(mu_lims), " target: ", n_county_final)
    return np.mean(mu_lims)

# def generate_forecasts_county(forecasting_parameters, total_ave, years):
#         # H_discrete,N,total,tau,years):
#     """
#     Generates the transition matrix for each year
#         :param forecasting_parameters:
#         :param df_blockgroup_county:
#         :param years:
#         :return:
#         """
#
#     # TODO: Switch to my hbar functions
#     # TODO: Write in function to calculate mu quickly
#
#     def forecast_county_single_race(H_discrete, total, taus, mu):
#
#         H_discrete_mu = H_discrete + mu * np.linspace(0, 1, len(H_discrete))
#         nonlocal years
#         if total>0 and tau>0 and len(H_discrete) > 1:
#             total = int(round(total))
#             transition_matrix = build_transition_matrix(H_discrete_mu,total, option = "right_multiply")
#             transition_matrix_per_interval = matrix_power(transition_matrix,int(round(tau*total*(years[1]-years[0]))))
#             forecast_county = np.zeros((len(years),total+1,total+1))
#             forecast_county[0,:,:] = transition_matrix
#             for i in range(1,len(years)):
#                 transition_matrix = np.matmul(transition_matrix_per_interval, transition_matrix)
#                 forecast_county[i, :, :] = transition_matrix
#         else:
#             forecast_county = np.nan*np.zeros((len(years),total+1,total+1))
#         return forecast_county
#     forecasts_county = {#"black": forecast_county_single_race(forecasting_parameters[forecasting_parameters["race"]=="black"]["H"],
#     #                        round(data_selected_blockgroup['black_00']),
#     #                        round(data_selected_blockgroup['total_00']),
#     #                        forecasting_parameters[forecasting_parameters["race"]=="black"]["H"]),
#                  # "hispanic": forecast_county_single_race(forecasting_parameters[forecasting_parameters["race"] == "hispanic"]["H"],
#                  #           round(data_selected_blockgroup['hispanic_00']),
#                  #           round(data_selected_blockgroup['total_00']),
#                  #           forecasting_parameters[forecasting_parameters["race"] == "hispanic"]["tau"]),
#                  "white": forecast_county_single_race(forecasting_parameters[forecasting_parameters["race"] == "white"]["H"].iloc[0],
#                                                       total_ave,
#                                                       forecasting_parameters[forecasting_parameters["race"] == "white"]["tau"].iloc[0],
#                                                       forecasting_parameters[forecasting_parameters["race"] == "white"]["mu"].iloc[0])
#                  }
#
#     return forecasts_county


def generate_blockgroup_forecast(race, mu, taus, hbars, data_selected_blockgroup, years, initial_year=2000):
    # H_discrete,N,total,tau,years):
    # Assumes that years is linearly increasing
    # TODO: Take care of when there is no matching forecast for that county
    def forecast_single_race(hbar, ns, N, total, mu, tau):
        H_discrete = interpolate_H_discrete_from_hbar(hbar, ns, total)
        H_discrete_mu = H_discrete + mu * np.linspace(0, 1, total + 1)
        nonlocal years
        if total > 0 and tau > 0 and len(H_discrete) > 1:
            total = round(total)
            transition_matrix = build_transition_matrix(H_discrete_mu, total, option="right_multiply")
            state_vector = np.zeros(total + 1)
            state_vector[N] = 1
            transition_matrix_per_interval = matrix_power(transition_matrix,
                                                          int(round(tau * total * (years[1] - years[0]))))
            forecast = np.zeros((len(years), total + 1))
            forecast[0, :] = state_vector
            for i in range(1, len(years)):
                state_vector = np.matmul(transition_matrix_per_interval, state_vector)
                forecast[i, :] = state_vector
        else:
            forecast = np.nan * np.zeros((len(years), total + 1))
        return forecast
    if isinstance(mu, pd.DataFrame):
        mu = mu[mu["race"] == race]["mu"].iloc[0]
    if initial_year == 2000:
        forecasts = forecast_single_race(hbars[hbars["race"] == race]["hbar"].iloc[0],
                                                   hbars[hbars["race"] == race]["ns"].iloc[0],
                                                   int(round(data_selected_blockgroup[race+'_00'].iloc[0])),
                                                   int(round(data_selected_blockgroup['total_00'].iloc[0])),
                                                   mu,
                                                   taus[taus["race"] == race]["tau"].iloc[0])
    elif initial_year == 2010:
        forecasts = forecast_single_race(hbars[hbars["race"] == race]["hbar"].iloc[0],
                                                   hbars[hbars["race"] == race]["ns"].iloc[0],
                                                   int(round(data_selected_blockgroup[race+'_10'].iloc[0])),
                                                   int(round(data_selected_blockgroup['total_10'].iloc[0])),
                                                   mu,
                                                   taus[taus["race"] == race]["tau"].iloc[0])

    return forecasts


# def generate_forecasts_selected_blockgroup(mus,tau,hbars, data_selected_blockgroup, years):
#         # H_discrete,N,total,tau,years):
#     # Assumes that years is linearly increasing
#     # TODO: Take care of when there is no matching forecast for that county
#     def forecast_single_race(hbar, ns, N, total, mu, tau):
#         H_discrete = interpolate_H_discrete_from_hbar(hbar, ns, total)
#
#         H_discrete_mu = H_discrete + mu * np.linspace(0, 1, total+1)
#         nonlocal years
#         if total>0 and tau>0 and len(H_discrete) > 1:
#             total = round(total)
#             transition_matrix = build_transition_matrix(H_discrete_mu,total, option = "right_multiply")
#             state_vector = np.zeros(total+1)
#             state_vector[N] = 1
#             transition_matrix_per_interval = matrix_power(transition_matrix,int(round(tau*total*(years[1]-years[0]))))
#             forecast = np.zeros((len(years),total+1))
#             forecast[0,:] = state_vector
#             for i in range(1,len(years)):
#                 state_vector = np.matmul(transition_matrix_per_interval, state_vector)
#                 forecast[i, :] = state_vector
#         else:
#             forecast = np.nan*np.zeros((len(years),total+1))
#         return forecast
#
#     forecasts = {"black": forecast_single_race(hbars[hbars["race"]=="black"]["hbar"].iloc[0],
#                             hbars[hbars["race"]=="black"]["ns"].iloc[0],
#                            int(round(data_selected_blockgroup['black_00'].iloc[0])),
#                            int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                           mus[mus["race"] == "black"]["mu"].iloc[0],
#                           tau_dict["black"], ),
#                  "hispanic": forecast_single_race(hbars[hbars["race"] == "hispanic"]["hbar"].iloc[0],
#                            hbars[hbars["race"] == "hispanic"]["ns"].iloc[0],
#                            int(round(data_selected_blockgroup['hispanic_00'].iloc[0])),
#                            int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                            mus[mus["race"] == "hispanic"]["mu"].iloc[0],
#                            tau_dict["hispanic"], ),
#                  "white": forecast_single_race(hbars[hbars["race"] == "white"]["hbar"].iloc[0],
#                         hbars[hbars["race"] == "white"]["ns"].iloc[0],
#                        int(round(data_selected_blockgroup['white_00'].iloc[0])),
#                        int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                        mus[mus["race"] == "white"]["mu"].iloc[0],
#                        tau_dict["white"], )}
#     return forecasts

# def generate_forecasts_selected_blockgroup_ARCHIVE(mus, tau_dict, hbars, data_selected_blockgroup, years):
#
# # H_discrete,N,total,tau,years):
# # Assumes that years is linearly increasing
# # TODO: Take care of when there is no matching forecast for that county
# def forecast_single_race(hbar, ns, N, total, mu, tau):
#     H_discrete = interpolate_H_discrete_from_hbar(hbar, ns, total)
#
#     H_discrete_mu = H_discrete + mu * np.linspace(0, 1, total + 1)
#     nonlocal years
#     if total > 0 and tau > 0 and len(H_discrete) > 1:
#         total = round(total)
#         transition_matrix = build_transition_matrix(H_discrete_mu, total, option="right_multiply")
#         state_vector = np.zeros(total + 1)
#         state_vector[N] = 1
#         transition_matrix_per_interval = matrix_power(transition_matrix,
#                                                       int(round(tau * total * (years[1] - years[0]))))
#         forecast = np.zeros((len(years), total + 1))
#         forecast[0, :] = state_vector
#         for i in range(1, len(years)):
#             state_vector = np.matmul(transition_matrix_per_interval, state_vector)
#             forecast[i, :] = state_vector
#     else:
#         forecast = np.nan * np.zeros((len(years), total + 1))
#     return forecast
#
# forecasts = {"black": forecast_single_race(hbars[hbars["race"] == "black"]["hbar"].iloc[0],
#                                            hbars[hbars["race"] == "black"]["ns"].iloc[0],
#                                            int(round(data_selected_blockgroup['black_00'].iloc[0])),
#                                            int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                                            mus[mus["race"] == "black"]["mu"].iloc[0],
#                                            tau_dict["black"], ),
#              "hispanic": forecast_single_race(hbars[hbars["race"] == "hispanic"]["hbar"].iloc[0],
#                                               hbars[hbars["race"] == "hispanic"]["ns"].iloc[0],
#                                               int(round(data_selected_blockgroup['hispanic_00'].iloc[0])),
#                                               int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                                               mus[mus["race"] == "hispanic"]["mu"].iloc[0],
#                                               tau_dict["hispanic"], ),
#              "white": forecast_single_race(hbars[hbars["race"] == "white"]["hbar"].iloc[0],
#                                            hbars[hbars["race"] == "white"]["ns"].iloc[0],
#                                            int(round(data_selected_blockgroup['white_00'].iloc[0])),
#                                            int(round(data_selected_blockgroup['total_00'].iloc[0])),
#                                            mus[mus["race"] == "white"]["mu"].iloc[0],
#                                            tau_dict["white"], )}
# return forecasts

def generate_forecasts_selected_blockgroup_multidecade(mus, tau_dict, hbars, data_selected_blockgroup, years):
    # H_discrete,N,total,tau,years):
    # Assumes that years is linearly increasing
    # TODO: Take care of when there is no matching forecast for that county
    def forecast_single_race(hbar, ns, N, total, tau, mu):
        H_discrete = interpolate_H_discrete_from_hbar(hbar, ns, total)
        H_discrete_mu = H_discrete + mu * np.linspace(0, 1, total + 1)
        nonlocal years
        if total > 0 and tau > 0 and len(H_discrete) > 1:
            total = round(total)
            state_vector = np.zeros(total + 1)
            state_vector[N] = 1
            years_initial=mus["year_initial"]
            years_final = mus["year_final"]
            mus_vector = mus["mu"]
            forecast = np.zeros((len(years), total + 1))
            forecast[0, :] = state_vector
            for (year_i, year_f) in zip(years_initial, years_final, mus_vector):
                H_discrete_mu = H_discrete + mu * np.linspace(0, 1, total + 1)
                transition_matrix = build_transition_matrix(H_discrete_mu, total, option="right_multiply")
                transition_matrix_per_interval = matrix_power(transition_matrix,
                                                              int(round(tau * total * (years[1] - years[0]))))
                years_mu = years[(years>year_i) & (years<=year_f)]
                for i in range(0, len(years_mu)):
                    state_vector = np.matmul(transition_matrix_per_interval, state_vector)
                    forecast[i, :] = state_vector
        else:
            forecast = np.nan * np.zeros((len(years), total + 1))
        return forecast
    print(data_selected_blockgroup)
    forecasts = {"black": forecast_single_race(hbars[hbars["race"] == "black"]["hbar"],
                                           hbars[hbars["race"] == "black"]["ns"],
                                           int(round(data_selected_blockgroup['black_00'])),
                                           int(round(data_selected_blockgroup['total_00'])),
                                           mus[mus["race"] == "black"]["mu"],
                                           tau_dict["black"], ),
             "hispanic": forecast_single_race(hbars[hbars["race"] == "hispanic"]["hbar"],
                                              hbars[hbars["race"] == "hispanic"]["ns"],
                                              int(round(data_selected_blockgroup['hispanic_00'])),
                                              int(round(data_selected_blockgroup['total_00'])),
                                              mus[mus["race"] == "hispanic"]["mu"],
                                              tau_dict["hispanic"], ),
             "white": forecast_single_race(hbars[hbars["race"] == "white"]["hbar"],
                                           hbars[hbars["race"] == "white"]["ns"],
                                           int(round(data_selected_blockgroup['white_00'])),
                                           int(round(data_selected_blockgroup['total_00'])),
                                           mus[mus["race"] == "white"]["mu"],
                                           tau_dict["white"], )}

    return forecasts


def interpolate_H_discrete_from_hbar(hbar,ns,total):
    hbar_func = interp1d(ns, hbar,kind='quadratic')
    ns = np.linspace(1/total,1-1/total,total-1)
    H_discrete = np.hstack((hbar_func(0),hbar_func(ns)-ns*np.log(ns) - (1-ns)*np.log(1-ns),hbar_func(1)))
    return H_discrete


def build_transition_matrix(H_discrete, total, option = "right_multiply"):
    n_s = np.linspace(0,1,total+1)
    dH_dN_increasing = total * (H_discrete[1:] - H_discrete[0:total])
    dH_dN_decreasing = total * (H_discrete[0:total] - H_discrete[1:])
    decreasing_N_rate = n_s[1:] / (1.0 + np.exp(dH_dN_decreasing))
    increasing_N_rate = (1.0 - n_s[0:total])/(1.0 + np.exp(dH_dN_increasing))
    constant_N_rate = 1.0 - np.append(0,decreasing_N_rate) - np.append(increasing_N_rate,0)
    if option == "left_multiply":
        transition_matrix = np.diagflat(constant_N_rate) + \
                            np.diagflat(decreasing_N_rate,-1) + \
                            np.diagflat(increasing_N_rate,1)
    elif option == "right_multiply":
        transition_matrix = np.diagflat(constant_N_rate) + \
                            np.diagflat(decreasing_N_rate,1) + \
                            np.diagflat(increasing_N_rate,-1)
    return transition_matrix


# def optimize_mu(N, total, H, max_steps, fraction_change,
#                   mu_min=None,
#                   mu_max=None,
#                   number_mus=20
#                   ):
#     # TODO: figure out a good guess for the maximum and minimum mu's
#     # TODO: Consider if I should optimize mu better by doing some sort of Gaussian thing instead of a simple list
#     mu_s = np.linspace(mu_min,mu_max,number_mus)
#     for mu in mu_s
#     return mu, tau

# def calculate_delta_V_steady_state(H, N, total, delta_n):
#     """
#     For quick forecasts, I want to find the value of delta_V for that decreases the expected
#     fraction of people across the county by delta_n.
#     :param H:
#     :param N:
#     :param total:
#     :param delta_n:
#     :return:
#     """
#     n_initial = sum(N)/sum(total)
#
#     return delta_V

#
# def calculate_delta_V_delta_t(H, N, total, delta_n, num_steps):
#     """
#     Calculates the correct value of delta_V so that the composition
#     of the county decreases by delta_n by the number of steps provided
#     by num_steps.
#     This is used to set a bound on the largest absolute value of delta_V so that
#     any larger would be unrealistic since it would decrease the composition too quickly.
#     :param H:
#     :param N:
#     :param total:
#     :param delta_n:
#     :param num_steps:
#     :return:
#     """
#     return delta_V
#
#
