# Overview

To extract Density-functional Fluctuation Theory (DFFT) functions quantifying the likelihood of observing racial/ethnic compositions in neighborhoods across the US, we employ a few computational strategies as described in an upcoming work. Code employing these strategies to both extract these functions and also use them to forecast compositional population dynamics into the future are gathered here.

## Computer Languages
All code here is written in either Python or Julia. Both are free, open-source languages. Heavy computational tasks are typically done in Julia while data processing and visualization are done in Python often using pandas, geopandas, and ipyleaflet modules.

## Data architecture
### Raw and processed data
All US census data were obtained through NHGIS.

1. Decadal data provide the exact population counts of persons by race and ethnicity within each block group in the US. Raw data is available in the raw-data folder. We further group this data by Asian NH (Non-Hispanic), Black NH, White NH, Hispanic, and other and provide that in the processed-data folder.

2. Longitudinal data were provided by NHGIS. NHGIS uses block level population counts and interpolation to approximate population counts within the same geometries demarcated in 2010 for both the 1990 and 2000 censuses. For more information on this process, see the `raw_data/nhgis0010_ts_geog2010_blck_grp_codebook.txt` file. We similarly group persons by race and ethnicity as described above and provide this data in the processed-data folder.

3. Traditional population projections used from https://doi.org/10.1038/sdata.2019.5 are used to forecast how neighborhoods are likely to change their compositions into the future. The referenced work provides five different forecasts, SSP1-SSP5. This data is processed by summing total populations across age cohorts and genders and saved in the processed-data folder.

## DFFT parameters
Code provided here take the population counts described above and produce the following DFFT parameters.
1. Global frustration functions, f(g). These functions are built by fitting discrete parameters for irregularly spaced compositions and then interpolating using a quadratic 1D spline fit. Then, to build the equations as represented in the paper, we add the function -(n)log(n) - (1-n)log(1-n) as described in the computational methods section of the paper. We provide here both the raw parameters as extracted from our functions and the interpolated version for a neighborhood size, s, of 1000 persons. Frustration functions can be found for Black/non-Black, Hispanic/non-Hispanic, and White/non-White for all three decades, 1990, 2000, 2010.

2. Vexation parameters, V^C. Alongside the optimization for the global frustration, we extract the Vexation for each county, each decade, and for each race/ethnicity (Black, Hispanic, White)

3. Time constants, τ. To find the time constant that optimizes the log-likelihood of the forecast, we simply generate forecasts (see below) for every census block group from the longitudinal population counts dataset for a range of values for the time constant τ. We do so using 1990 data/DFFT parameters to forecast 2000 data and also using 2000 data/DFFT parameters to forecast 2010 data. The optimal values are stored in `DFFT-parameters/taus_1990_2000.csv` and `DFFT-parameters/taus_2000_2010.csv`. We also provide the resulting forecasts using to calculate these optimum values at the OSF data repository associated with this project.

## Code for inferring DFFT functions
We use two methods for inferring the DFFT functions from population counts. This is a non-trivial, non-linear problem due to the different total population of each census block group and the underlying constraint that all counties possess the same segregation function, f(n), but are free to have their own values for Vexation. In summary, the first method, written in python, is faster but makes some mild assumptions to simplify the dataset. The second method, written in Julia, runs the exact optimization as described in our analysis paper, but is slow and works best when using an initial guess from the first method.

1. By binning the block groups into density domains and accounting for the appropriate combinatorial factors, we can more quickly approximate the optimum global frustration function and Vexation constants for each county. This function is found in the `DFFT-inference/BOptim.py` file alongside a Jupyter notebook `DFFT-inference/DFFT_fit_example.ipynb`. We find these functions to return parameters be very close to the best fit we obtain and thus is the best starting point to investigate the nature of these segregation functions.

2. To optimize DFFT functions further, we optimize the exact log-likelihood as presented in the paper. We assign discrete parameters for irregularly spaced compositions and use a gradient search method to maximum those results. This optimization is followed until the total log-likelihood across all block groups stops increasing past a certain threshold.

## Code for generating forecasts of compositional dynamics
To forecast how a neighborhood is likely to change it's composition into the future we perform the following step.
1. First, we must find the correct value for μ that shifts the composition of the county to match either the change in county composition observed for historical data, or the change projected for future forecasts. To do so, we build a state vector that matches the starting block group compositions, then use a binary search algorithm to find the correct value for μ that shifts the mean of that state vector the appropriate amount given the number of steps defined by the value assigned to the time constant τ.

2. Second, we iterate each neighborhood forward using the following metacode
    For each county
        For each race
            For each block group
                Start a new state vector representing the initial composition of that block group
                For each decade into the future
                    Calculate the new transition matrix using the value for μ
                    Iterate the state vector forward
                    Calculate the moments of the state vector probability distribution
                    If the forecast is of historical data, calculate the log-probability of observing the future distribution
                    Record the results
    This code is available in the two scripts `DFFT-forecast/forecast_historical.jl` for data used to validate forecasts or optimize the time scale τ and `DFFT-forecast/forecast_future.jl` for forecasts into the future using population projections. The results of these forecasts are saved into the folder `forecast-results`
