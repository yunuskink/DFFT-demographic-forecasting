#=
This script takes in a chosen population projection (here we use the 5 SSPs
from https://doi.org/10.1038/sdata.2019.5) to find the probability vector of
observing any given composition in each census block group across the US for the
years 2020-2100 starting from their block group compositions in 2010. This is
much slower than the simple 2020 forecast file due to the additional forecasts.

PART 1
Find the mu shift for each interval across all 5 SSPs and for each race. Save this as a dataframe

PART 2
For each state
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

## Part 1 #############################

##### LOAD IN DATA ######
counts_long_df = load_counts_longitudinal_df()
counts_long_df = counts_long_df[counts_long_df[:,:total_10].>0,:]
@load "./f_dict_V_df_hbar274.jld2" f_dict V_df ns
hbar_df = save_headaches_csv_hbar_uneven(f_dict,V_df,ns)
hbar_df = hbar_df[hbar_df[!,:decade].==2010,:]
@load "./taus_df.jld2" taus_df
projections_df = CSV.File("./data/population_projections_SSP2.csv") |> DataFrame
projections_df[!,[:year,:county_geoid]] = convert.(Int32,projections_df[!,[:year,:county_geoid]])
projections_df[:,:total] = sum(eachcol(projections_df[:,["1","2","3","4"]]))
rename!(projections_df,["1"=>"white","2"=>"black","3"=>"hispanic","4"=>"other"])
# Get only optimal taus obtained from 2000-2010 dynamics
# return taus_best_df
taus_best_df = get_best_taus(taus_df;decade_initial = 1990,decade_final = 2000, decade_H = 1990)

## Calculate mus
function generate_mu_from_projections()
    races = unique(hbar_df[!,:race])
    counties = unique(hbar_df[!,:county])
    counties = counties[1:5]
    years_initial = 2010:10:2090
    years_final = 2020:10:2100
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
                # println(row)
                append!(mus_df,DataFrame(mu_row))
            end
        end
        println("County: ", county)
    end
    return mus_df
end
# using Profile


mus_df=generate_mu_from_projections()
@save "mus_df_27_01_2021.jld2" mus_df
## Forecast blockgroups
# Ideas to store reults
# forecasted_mean_df: Just keeps track of the means/variances/skewnesses. Each decade/race combo
# should have its own column, making 27*3 columns. I'm fine with that. Populate it first with Missings
# and then fill them in as things go.
# forecasted_probability_df: One file for each state/race combo. This should keep files under 1 GB. Initialize within loop
# Initialize dataframe: Only initialize forecasted moments once
year_strings = [string(d) for d in 2020:10:2100]
moment_strings = ["mean", "variance", "skewness"]
colnames = [race*"_"*year*"_"*moment for race in races for year in year_strings for moment in moment_strings]
forecasted_moments_df = DataFrame(geoid_blockgroup = counts_long_df[:,:geoid_blockgroup])
for colname in colnames
    forecasted_moments_df[:,colname] .= NaN64
end

function generate_forecasts_2100()
    counties_all = unique(hbar_df[!,:county])
    counties_all = unique(mus_df[!,:geoid_county])
    counties_all = counties_all[2]
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
                        # probability_df_row[!,string(year_f)].=state_vector
                        forecasted_probability_df[forecasted_probability_df[:,:geoid_blockgroup].==blockgroup,string(year_f)] = [state_vector]
                    end
                    # append!(forecasted_probability_df,probability_df_row)
                end
                println("County: ", county)
            end
            @save "./"*race*"/"*string(1)*"_forecasted_probability_df_27_01_2021.jld2" forecasted_probability_df
        end
        @save "forecasted_moments_df_27_01_2021.jld2" forecasted_moments_df
    end
end

@time generate_forecasts_2100()
