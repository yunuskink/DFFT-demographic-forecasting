## Apps overview
In order to make forecasts available to the public, we publish interactive apps here that allow the user to explore forecasts that DFFT generates. Specifically, we provide two separate apps, one to check the accuracy of forecasts of historic data and another to explore forecasts of unreleased data. These apps are written in python and can be run in the cloud without needing to download or install software, however it will typically increase the speed if users download code and run it on their machine locally.

### Instructions
1. Use this link to open the notebook, or download this repository and run the `validation_and_forecasting_apps.ipynb`.
2. Click `Run all` from the `Cell` dropdown menu to run all the cells in the notebook.


### Historic validation app
To explore how well DFFT forecasts already released data, we present the results of the forecast of 2010 census block group compositions using 1990 and 2000 data. To do so...
3. Use the search tool or drag the marker to a location of interest. User will then be presented with a choropleth of the county in which the marker is placed. Choose the desired layer to view in the choropleth.
4. On the right-hand side, user can drag the year slider at the bottom to see how the probability of observing a given composition in the neighborhood is expected to evolve into the future. When set to the year 2000, this probability distribution has a value of 1 at the known composition observed in the year 2000. As time goes on, this distribution will spread. When the slider is set to the year 2010, the user can then compare the forecasted probability distribution to the observed composition in the year 2010.

### Future forecast app
For forecasting future compositions, we do not have access to exactly how county compositions will change. Instead, we must resort to traditional population projections of coarse county composition changes in order to use DFFT to forecast how individual block groups are likely to change into the future.

1. Use this link to open the notebook, or download this repository and run the `future_forecast_GUI.ipynb`.
2. Click `Run all` from the `Cell` dropdown menu to run all the cells in the notebook.
