#=
This script takes in a chosen population projection (here we use the 5 SSPs
from https://doi.org/10.1038/sdata.2019.5) to find the probability vector of
observing any given composition in each census block group across the US for the
year 2020 starting from their block group compositions in 2010.

PART 1
Find the mu shift for the chosen population projection (SSP) for each race. Save this as a dataframe

PART 2
For each state
    For each county
        For each race
            For each blockgroup
                Start a new state vector
                For each interval
                    Calculate the new transition matrix
                    iterate the state vector forward
                    Add the state vector to the dataframe

STRUCTURE OF RESULTS
mu_df:
Columns -> year_initial
           composition_initial
           year_finial
           composition_final
           mu

forecasted_probability_df: This is currently not generated since the size is so
large, but simply uncommenting labelled lines below will be able to produce these
dataframes for each state to limit file sizes.
Columns -> geoid_blockgroup
           year
           probability

forecasted_moments_df:
This is a more manageable dataframe that keeps only the mean, variance, and skewness of the forecasts.
=#
include("generate_forecasts_subfunctions.jl")
include("file_import_export.jl")

## Part 1 #############################
SSP = 2
projections_df = CSV.File("./data/population_projections_SSP"*string(SSP)*".csv") |> DataFrame
@load "./taus_df.jld2" taus_df
taus_best_df = get_best_taus(taus_df;decade_initial = 2000,decade_final = 2010, decade_H = 2000)

years_initial = 2010
years_final = 2020

##### LOAD IN DATA ######
counts_long_df = load_counts_longitudinal_df()
counts_long_df = counts_long_df[counts_long_df[:,:total_10].>0,:]
@load "./f_dict_V_df_hbar274.jld2" f_dict V_df ns
hbar_df = save_headaches_csv_hbar_uneven(f_dict,V_df,ns)
hbar_df = hbar_df[hbar_df[!,:decade].==2010,:]
projections_df[!,[:year,:county_geoid]] = convert.(Int32,projections_df[!,[:year,:county_geoid]])
projections_df[:,:total] = sum(eachcol(projections_df[:,["1","2","3","4"]]))
rename!(projections_df,["1"=>"white","2"=>"black","3"=>"hispanic","4"=>"other"])

## Calculate mus
function generate_mu_from_projections()
    races = unique(hbar_df[!,:race])
    counties = unique(hbar_df[!,:county])
    counties = counties[1:5]
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
        projections_subset = projections_df[projections_df[!,:county_geoid].==county,:]
        totals = counts_subset[:,"total_10"]
        total_mean = convert(Int32,round(sum(totals)/length(totals)))
        for race in races
            # Start a state vector from 2010 data
            Ns = counts_subset[:,race*"_10"]
            state_vector = calculate_state_vector_initial(Ns,totals,total_mean)
            # Get the discrete form of H
            H = calculate_H_discrete(hbar_subset[(hbar_subset[!,:race].==race),[:hbar,:ns]],total_mean)
            tau = taus_subset[taus_subset[!,:race].==race,:tau][1]
            for (year_i,year_f) in zip(years_initial,years_final)
                fraction_final = projections_subset[projections_subset[!,:year].==year_f,race]/
                                 projections_subset[projections_subset[!,:year].==year_f,:total]
                # println(projections_subset[projections_subset[!,:year].==year_f,:])
                mu,state_vector = calculate_mu(state_vector,fraction_final[1],tau,total_mean,H;mu_lims = [-2.0 2.0],)
                mu_row = DataFrame(race = race, geoid_county = county, tau = tau,
                                mu=mu, year_initial = year_i, year_final=year_f)
                append!(mus_df,DataFrame(mu_row))
            end
        end
        println("County: ", county)
    end
    return mus_df
end

mus_df=generate_mu_from_projections()
@save "mus_df_2020.jld2" mus_df
## Forecast blockgroups
# Ideas to store reults
races = unique(hbar_df[!,:race])
year_strings = [string(d) for d in 2020]
moment_strings = ["mean", "variance", "skewness"]
colnames = [race*"_"*year*"_"*moment for race in races for year in year_strings for moment in moment_strings]
forecasted_moments_df = DataFrame(geoid_blockgroup = counts_long_df[:,:geoid_blockgroup])
for colname in colnames
    forecasted_moments_df[:,colname] .= NaN64
end

function generate_forecasts_2020()
    counties_all = unique(hbar_df[!,:county])
    # counties_all = unique(mus_df[!,:geoid_county])
    counties_all = counties_all[1]
    states = convert.(Int32,unique(floor.(counties_all./1000)))
    for race in races
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
                projections_subset = projections_df[projections_df[!,:county_geoid].==county,race]
                (geoids_blockgroup,Ns,totals) = eachcol(counts_long_df[counts_long_df[:,:geoid_county].==county,["geoid_blockgroup",race*"_10","total_10"]])
                mus_subset = mus_df[(mus_df[:,:geoid_county].==county).&
                                    (mus_df[:,:race].==race),:]
                for (blockgroup,N,total) in zip(geoids_blockgroup,Ns,totals)
                    # probability_df_row = DataFrame(geoids_blockgroup=blockgroup)
                    # for year_string in year_strings
                    #     probability_df_row[:,year_string] .= Any[]
                    # end
                    state_vector = zeros(Float64,total+1)
                    state_vector[N+1] = 1
                    # Get the discrete form of H
                    H = calculate_H_discrete(hbar,total)
                    # println(state_vector)
                    # println(total)
                    num_absolute_steps = convert(Int64,round(total*tau*10))
                    for (year_i,year_f) in zip(years_initial,years_final)
                        mu = mus_subset[mus_subset[:,:year_final].==year_f,:mu][1]
                        H_mu = H.+mu.*(0:1/total:1)
                        transition_matrix = calculate_transition_matrix(H_mu,total)
                        for i=1:num_absolute_steps
                            state_vector = transition_matrix * state_vector
                        end
                        mn = sum((0:1/total:1).*state_vector)
                        var = sum(((0:1/total:1).^2).*state_vector) - mn^2
                        skew = sum((((0:1/total:1).-mn).^3).*state_vector)/(var^(3/2))
                        idx = findfirst(forecasted_moments_df[:,:geoid_blockgroup].==blockgroup)
                        forecasted_moments_df[idx,[race*"_"*string(year_f)*"_mean",
                                                    race*"_"*string(year_f)*"_variance",
                                                    race*"_"*string(year_f)*"_skewness"]]=(mn,var,skew)
                        forecasted_probability_df[forecasted_probability_df[:,:geoid_blockgroup].==blockgroup,string(year_f)] = [state_vector]
                    end
                    # Uncomment the following line to build the full probability distributions
                    # append!(forecasted_probability_df,probability_df_row)
                end
                println("County: ", county)
            end
            mkpath("./SSP"*string(SSP)*"/"*race*"/"*string(state))
            @save "./SSP"*string(SSP)*"/"*race*"/"*string(state)*"/forecasted_probability_2010_2020.jld2" forecasted_probability_df
        end
        @save "./forecasted_moments/forecast_2010_2020.jld2" forecasted_moments_df
    end
end

@time generate_forecasts_2020()

# state = 1
# @load "./"*race*"/"*string(state)*"/forecasted_probability_2010_2020.jld2" forecasted_probability_df

## Writing CSV files
using JLD2
using DataFrames
using CSV

@load "./forecasted_moments/2010_to_2020/SSP1/forecast_2010_2020.jld2" forecasted_moments_df
CSV.write("./forecasted_moments/2010_to_2020/SSP1/forecast_2010_2020.csv",forecasted_moments_df)
@load "./forecasted_moments/2010_to_2020/SSP2/forecast_2010_2020.jld2" forecasted_moments_df
CSV.write("./forecasted_moments/2010_to_2020/SSP2/forecast_2010_2020.csv",forecasted_moments_df)
@load "./forecasted_moments/2010_to_2020/SSP3/forecast_2010_2020.jld2" forecasted_moments_df
CSV.write("./forecasted_moments/2010_to_2020/SSP3/forecast_2010_2020.csv",forecasted_moments_df)
@load "./forecasted_moments/2010_to_2020/SSP4/forecast_2010_2020.jld2" forecasted_moments_df
CSV.write("./forecasted_moments/2010_to_2020/SSP4/forecast_2010_2020.csv",forecasted_moments_df)
@load "./forecasted_moments/2010_to_2020/SSP5/forecast_2010_2020.jld2" forecasted_moments_df
CSV.write("./forecasted_moments/2010_to_2020/SSP5/forecast_2010_2020.csv",forecasted_moments_df)
