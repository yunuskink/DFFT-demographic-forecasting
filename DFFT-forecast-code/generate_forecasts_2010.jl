#=
Going to try assigning a value for Tau. Really,

SHITSHIT! I fixed the goddamn so the R^2 for white is pretty damn good, but now Hispanic
is awful and Black is just OK. Something is fucked up about Hispanic at least. The
forecast for Hispanic

SHITSHITSHIT! R^2 for white people is negative between DFFT and absolute... It is
a little better than proportional. For Balck and Hispanic, however, everything seems
like what I expected... WTF!!!

Steps to try in a different script. Plot the residuals as a function of the initial composition
for all 3 (absolute, DFFT, proportional). See where the fuck this is failing.

generate_forecasts_2010 is different from generate_forecasts_2100 because, since we
know the resulting compositions, we can additionally calculate the log-likelihood
of the forecast. I will also calculate further statistics here for convenience later
Specifically, I'll get the residuals from DFFT, proportional, and unchanging.

PART 1
Find the mu shift for each interval across all 5 SSPs and for each race. Save this as a dataframe

PART 2
For each county
    For each race
        For each blockgroup
            Start a new state vector
            For each interval
                Calculate the new transition matrix
                iterate the statevector forward
                Add the statevector to the dataframe

STRUCTURE OF RESULTS
mu_df:
Columns -> year_initial
           composition_initial
           year_finial
           composition_final
           mu

forecast_df:
Columns -> geoid_blockgroup
           year
           probability

hbar_df:

=#
include("generate_forecasts_subfunctions.jl")
include("file_import_export.jl")
using Random
## Part 1 #############################

##### LOAD IN DATA ######
counts_long_df = load_counts_longitudinal_df()
counts_long_df = counts_long_df[(counts_long_df[:,:total_10].>0).&
                                (counts_long_df[:,:total_00].>0).&
                                (counts_long_df[:,:total_90].>0),:]
@load "./f_dict_V_df_hbar274.jld2" f_dict V_df ns
hbar_df = save_headaches_csv_hbar_uneven(f_dict,V_df,ns)
hbar_df = hbar_df[hbar_df[!,:decade].==2000,:]
@load "./taus_df.jld2" taus_df

# Get only optimal taus obtained from 2000-2010 dynamics
# return taus_best_df
taus_best_df = get_best_taus(taus_df;decade_initial = 1990,decade_final = 2000, decade_H = 1990)

## Calculate mus

function generate_mu_from_projections()
    races = unique(hbar_df[!,:race])
    counties = unique(hbar_df[!,:county])
    # counties = counties[1:5]
    years_initial = 2000
    years_final = 2010
    mus_df = DataFrame(race = String[],geoid_county = Int32[],
                        year_initial = Int32[],
                        year_final = Int32[],
                        tau = Float64[],
                        mu = Float64[])

    for county in counties
        # Subset all the dataframes
        taus_subset = taus_best_df[taus_best_df[!,:geoid_county].==county,:]
        counts_subset = counts_long_df[counts_long_df[!,:geoid_county].==county,:]
        hbar_subset = hbar_df[hbar_df[!,:county].==county,:]
        totals = counts_subset[:,"total_00"]
        total_mean = convert(Int32,round(sum(totals)/length(totals)))
        for race in races
            # Start a state vector from 2010 data
            Ns = counts_subset[:,race*"_00"]
            state_vector = calculate_state_vector_initial(Ns,totals,total_mean)
            # Get the discrete form of H
            H = calculate_H_discrete(hbar_subset[(hbar_subset[!,:race].==race),[:hbar,:ns]],total_mean)
            tau = taus_subset[taus_subset[!,:race].==race,:tau][1]
            for (year_i,year_f) in zip(years_initial,years_final)
                fraction_final = sum(counts_subset[:,race*"_10"])/sum(counts_subset[:,"total_10"])
                mu,state_vector = calculate_mu(state_vector,fraction_final[1],tau,total_mean,H;mu_lims = [-2.0 2.0],)
                mu_row = DataFrame(race = race, geoid_county = county, tau = tau,
                                mu=mu, year_initial = year_i, year_final=year_f)
                # println(row)
                append!(mus_df,DataFrame(mu_row))
            end
        end
        println("County: ", county)
    end
    return mus_df
end
# using Profile

# @time mus_df=generate_mu_from_projections()
# @save "mus_df_validation_27_01_2021.jld2" mus_df

## Forecast blockgroups
# Ideas to store reults
# forecasted_mean_df: Just keeps track of the means/variances/skewnesses. Each decade/race combo
# should have its own column, making 27*3 columns. I'm fine with that. Populate it first with Missings
# and then fill them in as things go.
# forecasted_probability_df: One file for each state/race combo. This should keep files under 1 GB. Initialize within loop
# Initialize dataframe: Only initialize forecasted moments once
@load "mus_df_validation_27_01_2021.jld2" mus_df
year_strings = [string(d) for d in 2010]
moment_strings = ["mean", "variance", "skewness", "lglik", "res_sq_DFFT","res_sq_proportional","res_sq_absolute"]
# races = unique(hbar_df[:,:race])

function generate_forecasts_2010()
    counties_all = unique(hbar_df[!,:county])
    # counties_all = counties_all[randperm(length(counties_all))[1:20]]
    # counties_all = unique(mus_df[!,:geoid_county])
    # counties_all = counties_all[1:10]
    states = convert.(Int32,unique(floor.(counties_all./1000)))
    year_f = 2010
    for race in races
        forecasted_moments_df = DataFrame(geoid_blockgroup = counts_long_df[:,:geoid_blockgroup])
        colnames = [race*"_"*year*"_"*moment for year in year_strings for moment in moment_strings]
        for colname in colnames
            forecasted_moments_df[:,colname] .= NaN64
        end
        for state in states
            # forecasted_probability_df = DataFrame(geoid_blockgroup = Int64[])
            forecasted_probability_df = DataFrame(geoid_blockgroup = counts_long_df[
                                                    floor.(counts_long_df[:,:geoid_county]./1000).==state,:geoid_blockgroup])
            for year_string in year_strings
                forecasted_probability_df[:,year_string] .= [Any[]]
            end
            counties = counties_all[floor.(counties_all./1000).==state]
            for county in counties
                tau = taus_best_df[(taus_best_df[!,:geoid_county].==county).&
                                            (taus_best_df[!,:race].==race),:tau][1]
                hbar = hbar_df[(hbar_df[!,:county].==county).&(hbar_df[!,:race].==race),[:hbar,:ns]]
                (geoids_blockgroup,Ns,totals,Ns_10,totals_10) = eachcol(counts_long_df[counts_long_df[:,:geoid_county].==county,
                                                        ["geoid_blockgroup",race*"_00","total_00",race*"_10","total_10",]])
                counts_subset = counts_long_df[(counts_long_df[:,:geoid_county].==county),:]
                mu = mus_df[(mus_df[:,:year_final].==year_f).&(mus_df[:,:geoid_county].==county).&(mus_df[:,:race].==race),:mu][1]
                fraction_initial = sum(counts_subset[:,race*"_00"])/sum(counts_subset[:,"total_00"])
                fraction_final = sum(counts_subset[:,race*"_10"])/sum(counts_subset[:,"total_10"])
                if fraction_initial>fraction_final #Decreased in composition
                    Ns_prop = Ns.*(fraction_final/fraction_initial)
                else
                    Ns_prop = Ns .+ ((totals.-Ns).*((1-fraction_final)/(1-fraction_initial)))
                end
                for (blockgroup,N,total,N_10,total_10,N_prop) in zip(geoids_blockgroup,Ns,totals,Ns_10,totals_10,Ns_prop)
                    state_vector = zeros(Float64,total+1)
                    state_vector[N+1] = 1
                    num_absolute_steps = convert(Int64,round(total*tau*10))
                    H_mu = calculate_H_discrete(hbar,total,mu=mu)
                    transition_matrix = calculate_transition_matrix(H_mu,total)
                    for i=1:num_absolute_steps
                        state_vector = transition_matrix * state_vector
                    end
                    N_10_adj = convert(Int32,round(N_10*total/total_10))
                    mn = sum((0:1/total:1).*state_vector)
                    var = sum(((0:1/total:1).^2).*state_vector) - mn^2
                    skew = sum((((0:1/total:1).-mn).^3).*state_vector)/(var^(3/2))
                    lglik = log(state_vector[N_10_adj+1])
                    res_sq_DFFT = (total*mn - N_10_adj)^2
                    res_sq_absolute = (N - N_10_adj)^2
                    res_sq_proportional = (N_prop - N_10_adj)^2
                    idx = findfirst(forecasted_moments_df[:,:geoid_blockgroup].==blockgroup)
                    forecasted_moments_df[idx,[race*"_"*string(year_f)*"_mean",
                                                race*"_"*string(year_f)*"_variance",
                                                race*"_"*string(year_f)*"_skewness",
                                                race*"_"*string(year_f)*"_lglik",
                                                race*"_"*string(year_f)*"_res_sq_DFFT",
                                                race*"_"*string(year_f)*"_res_sq_absolute",
                                                race*"_"*string(year_f)*"_res_sq_proportional",
                                                ]]=(mn,var,skew,lglik,res_sq_DFFT,res_sq_absolute,res_sq_proportional)
                    # probability_df_row[!,string(year_f)].=state_vector
                    forecasted_probability_df[forecasted_probability_df[:,:geoid_blockgroup].==blockgroup,string(year_f)] = [state_vector]
                    # append!(forecasted_probability_df,probability_df_row)
                end
                println("County: ", county)
            end
            @save "/media/yunuskink/1 TB HDD/forecasts/"*race*"/"*string(state)*"_forecasted_probability_validation_df_27_01_2021.jld2" forecasted_probability_df
            @save "/media/yunuskink/1 TB HDD/forecasts/"*race*"/forecasted_moments_validation_df_27_01_2021.jld2" forecasted_moments_df
        end
    end
end


#
@time generate_forecasts_2010()
#
# # Combine the 3 races moments dataframes
@load "/media/yunuskink/1 TB HDD/forecasts/"*"white"*"/forecasted_moments_validation_df_27_01_2021.jld2" forecasted_moments_df
forecasted_moments_df_white = forecasted_moments_df
@load "/media/yunuskink/1 TB HDD/forecasts/"*"black"*"/forecasted_moments_validation_df_27_01_2021.jld2" forecasted_moments_df
forecasted_moments_df_black = forecasted_moments_df
@load "/media/yunuskink/1 TB HDD/forecasts/"*"hispanic"*"/forecasted_moments_validation_df_27_01_2021.jld2" forecasted_moments_df
forecasted_moments_df_hispanic = forecasted_moments_df
forecasted_moments_df = innerjoin(forecasted_moments_df_black,forecasted_moments_df_hispanic,on=:geoid_blockgroup)
forecasted_moments_df = innerjoin(forecasted_moments_df,forecasted_moments_df_white,on=:geoid_blockgroup)
# using CSV
CSV.write("/media/yunuskink/1 TB HDD/forecasts/forecasted_moments_validation_df_27_01_2021.csv", forecasted_moments_df)

## Adding on columns to forecasted_moments_df for visualization GUI
# if layer_type=='2000':
#     name = "fraction_" + race + "_" + layer_type
# elif layer_type == '2010':  # Choosing just a decade, not the dynamic change
#     name = "fraction_" + race + "_" + layer_type
# elif layer_type == '2000->2010 Observed':  # Choosing just a decade, not the dynamic change
#     name = "delta_fraction_observed_" + race + "_2000_2010"
# elif layer_type == '2000->2010 Forecasted':  #Choosing just a decade, not the dynamic change
#     name = "delta_fraction_forecasted_" + race + "_2000_2010"
# elif layer_type == '2000->2010 Error':  # Choosing just a decade, not the dynamic change
#     name = "error_" + race + "_2000_2010"
# elif layer_type == '2000->2010 Log-likelihood':  # Choosing just a decade, not the dynamic change
#     name = "error_" + race + "_2000_2010"

forecasted_moments_df = innerjoin(forecasted_moments_df,counts_long_df,on=:geoid_blockgroup)

races = ["black","hispanic","white"]
for race in races
    name = "fraction_" * race * "_" * "1990"
    forecasted_moments_df[:,name] = forecasted_moments_df[:,race*"_90"]./forecasted_moments_df[:,"total_90"]
    name = "fraction_" * race * "_" * "2000"
    forecasted_moments_df[:,name] = forecasted_moments_df[:,race*"_00"]./forecasted_moments_df[:,"total_00"]
    name = "fraction_" * race * "_" * "2010"
    forecasted_moments_df[:,name] = forecasted_moments_df[:,race*"_10"]./forecasted_moments_df[:,"total_10"]
    name = "delta_fraction_observed_" * race * "_2000_2010"
    forecasted_moments_df[:,name] = forecasted_moments_df[:,"fraction_" * race * "_" * "2010"].-forecasted_moments_df[:,"fraction_" * race * "_" * "2000"]
    name = "delta_fraction_forecasted_" * race * "_2000_2010"
    forecasted_moments_df[:,name] = forecasted_moments_df[:,race * "_2010_mean"].-forecasted_moments_df[:,"fraction_" * race * "_" * "2000"]
    name = "error_" * race * "_2000_2010"
    forecasted_moments_df[:,name] = (forecasted_moments_df[:,race * "_2010_mean"].-forecasted_moments_df[:,"fraction_" * race * "_" * "2010"])./forecasted_moments_df[:,race * "_2010_variance"]
    name = "lglik_" * race * "_2000_2010"
    forecasted_moments_df[:,name] = (forecasted_moments_df[:,race * "_2010_lglik"])./forecasted_moments_df[:,"total_00"]
end
using CSV
5+5
CSV.write("/media/yunuskink/1 TB HDD/forecasts/forecasted_moments_validation_df_27_01_2021.csv", forecasted_moments_df)

histogram(forecasted_moments_df[:,:error_white_2000_2010])

CSV.write("/media/yunuskink/1 TB HDD/forecasts/parameters/hbar_df.csv", hbar_df)
CSV.write("/media/yunuskink/1 TB HDD/forecasts/parameters/taus_df_validation.csv", taus_best_df)
CSV.write("/media/yunuskink/1 TB HDD/forecasts/parameters/mus_df_validation.csv", mus_df)
