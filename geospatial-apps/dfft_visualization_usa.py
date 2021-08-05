import ipyleaflet
from branca.colormap import linear
from shapely.geometry import shape, Point
import numpy as np
import pandas as pd
import ipywidgets as widgets
import bqplot as bq
import dfft_forecasting as dfft_fcast
import dfft_io

# Just do everything for decades, not every 5 years

"""
Here are functions for an interactive app that allows users to find a neighborhood in the US using 
addresses and forecast how that neighborhood will change its composition from the year 2010 onwards
given a population projection of how a given county changes its composition over that same time period.

This app is built in python using ipywidgets and ipylealet and the app is run using an accompanying
Jupyter notebook   
"""

# TODO: Automatically uncheck "Custom county forecast" checkbox when county changes
# TODO: Add in a year option for the custom county forecast
# TODO: Show the county population projection somewhere? Probably need a button to switch between the two
# TODO: Add in an option to switch between the 5 SSPs
# TODO: Add a footer that outputs important notes, like the total population of the neighborhood, mu, etc...
# TODO: Figure out the bug causing the choropleth to be mismatched with the frame

# TODO: Create another GUI that does the full county probability heatmap
# TODO: Redo old code as a validation app. Keep the two separate for now.
# TODO: Update choropleths to see instead the forecasts of the mean change.


# INPUTS NEEDED
# TODO: build mus_df. each county/race combo should have it's own row
# TODO: build forecasted_mean_df. Each block group/race has its own row. Then a certain column is chosen for choropleth
# TODO: build taus_00_10_df, just selecting the best tau for each county and each race
# TODO: build hbars_10_df and hbars_00_df, just saving it as a csv from Julia


def update_county_data(df_blockgroup, geoid_county):
    df_blockgroup_county = df_blockgroup[df_blockgroup['geoid_county']==geoid_county]
    df_blockgroup_county.replace([np.inf, -np.inf], np.nan)
    df_blockgroup_county = df_blockgroup_county.dropna(inplace=False)
    geo_df_blockgroup, geojson_blockgroup = dfft_io.load_blockgroup_GIS(geoid_county)
    # Remove from shapefiles any block groups that are not present in df_blockgroup_county
    rows_feature = [feature['properties']['GEOID10'] for feature in geojson_blockgroup['features']]
    rows_df_blockgroup_county = [str(int(row['geoid_blockgroup'])) for index, row in df_blockgroup_county.iterrows()]
    blockgroups_to_remove = np.setxor1d(rows_feature, rows_df_blockgroup_county, assume_unique=True)
    indices_to_remove = []
    for i in range(len(geojson_blockgroup['features'])):
        if geojson_blockgroup['features'][i]['id'] in blockgroups_to_remove:
            # if geojson_blockgroup['features'][i]['GEOID10'] in blockgroups_to_remove:
            indices_to_remove.append(i)
    for index in sorted(indices_to_remove, reverse=True):
        geojson_blockgroup['features'].pop(index)
    blockgroups_intersect = np.intersect1d(geo_df_blockgroup['GEOID10'],df_blockgroup_county['geoid_blockgroup'])
    # print(blockgroups_intersect)
    return geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup



def build_forecasting_app_old(forecasting_df_filename = "./data/forecasting_parameters_df_Boris.csv"):
    """
    Users will be able to see how the mean composition of a county will change over time
    given their choice of a population projection. They can also define a county composition shift
    and see the same results except they won't see the way the entire county changes since that is
    too slow.

    First, user will need to use the search tool to place the marker within a county of interest.
    Second, user will be presented with an animated heatmap for the probability of observing a
    future composition given the initial composition for a neighborhood of 1000 persons
    Third, user can choose from a dropdown menu to see 9 separate entries for the 2000,2010 observed
    and 2010 forecasted mean for the selected county.
    Fourth, user will be able to click on a neighborhood of interest and see an animated plot for
    the dynamics of that specific neighborhood
    To explore a new county, the user need only use the search tool to find a new county.
    BACKEND:
    First, I need a function that is called when the search tool is used. A county should then be
    selected and the shapefile, Headache, optimum delta_V, and optimum delta_T should be loaded
    from file or calculated. This will be a dataframe.
    Second,

    App widgets...
    1) race_selector: lets user choose whether they want Black/non-Black,White/non-White,Hispanic/non-Hispanic
    2) decade_selector: lets user choose whteher to see just a single decade or the change between two decades
    3) SearchControl: Lets users find new location
    4) marker: Used to collect the location of a given county for analysis
    5) year_slider: Lets users forecast further into the future and intervening years

    :return:
    """
    initial_position = (41.7406, -87.7093) #Somewhere in Cook county
    years = np.linspace(2020, 2100, 9, 'Int')

    ###### LOADING DATA ################
    df_blockgroup = pd.read_pickle("./data/blockgroup_longitudinal_data")
    df_blockgroup = dfft_io.add_columns_to_blockgroup_df(df_blockgroup) #Calculate columns like fractional change
    geo_df_county = dfft_io.load_county_GIS()
    hbars_df = dfft_io.load_dataframe_with_hbar_from_csv("./parameters/hbar_df.csv")
    hbars_df = hbars_df[hbars_df["decade"]==2010]# hbars for 2010
    mus_df = pd.read_csv("./parameters/mus_df_SSP2.csv")  # mus used to forecast 2020-2100. Also used for the forecasts
    taus_df = pd.read_csv("./parameters/taus_df_forecasting.csv")  # taus obtained from 2000-2010 dynamics
    df_blockgroup = pd.read_csv("./data/forecasted_moments_SSP2_df.csv")  # each block group mean 2020-2100
    ######### FORECAST FOR INITIAL NEIGHBORHOOD ################
    # TODO: Consider if I should bother with the initial neighborhood. Maybe it would be best to have the user do something first like move a marker or search
    # Get id for initial county
    geoid_county = geo_df_county.loc[geo_df_county.contains(Point(initial_position[1], initial_position[0]))]
    geoid_county = geoid_county.iloc[0]['GEOID']
    print(mus_df)
    # Get the parameters for that county
    mus = mus_df[mus_df["geoid_county"]==geoid_county]
    taus = taus_df[taus_df["geoid_county"]==geoid_county]
    hbars = hbars_df[hbars_df["geoid_county"]==geoid_county]

    # Find the blockgroup that matches the location
    geo_df_blockgroup, geojson_blockgroup = dfft_io.load_blockgroup_GIS(geoid_county)
    geoid_blockgroup = geo_df_blockgroup.loc[
        geo_df_blockgroup.contains(Point(initial_position[1], initial_position[0]))]
    geoid_blockgroup = geoid_blockgroup.iloc[0]['GEOID10']
    data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]

    # Generate forecasts of that neighborhood
    forecasts = dfft_fcast.generate_forecasts_selected_blockgroup_multidecade(mus, taus, hbars, data_selected_blockgroup, years)
    ######## WIDGETS ############
    # race_selector = widgets.Dropdown(options=('black', 'hispanic', 'white'),
    race_selector = widgets.ToggleButtons(options=(['black','hispanic','white']),
                                          description = 'Race/Ethnicity:',
                                        # layout=widgets.Layout(width='20%'),
                                          style={'description_width': 'initial'},
                                          disabled=True,
                                          )

    choropleth_selector = widgets.ToggleButtons(options=(['Raw mean compositions','Mean change in composition']),
                                          description = 'Choropleth option:',
                                        # layout=widgets.Layout(width='20%'),
                                          style={'description_width': 'initial'},
                                          disabled=True)
    decade_selector = widgets.Dropdown(options=('2010', '2020', '2030','2040','2050','2060','2070','2080','2090','2100'),
                                       description='Year:',
                                       # layout=widgets.Layout(width='20% '),
                                       style={'description_width': 'initial'})
    projection_selector = widgets.Dropdown(options=('SSP1', 'SSP2', 'SSP3', 'SSP4', 'SSP5'),
                                           description='County projection:',
                                           # layout=widgets.Layout(width='20% '),
                                           style={'description_width': 'initial'})
    custom_forecast_checkbox = widgets.Checkbox(
                                        value=False,
                                        description='Custom county forecast',
                                        disabled=True,
                                        indent=False)
    custom_forecast_text_initial = widgets.BoundedFloatText(
                                                value=0,
                                                min=-1.0,
                                                max=1.0,
                                                step=0.01,
                                                description='Change in county composition:',
                                                disabled=True,
                                                continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,
                                                readout_format='.1f',
                                                style={'description_width': 'initial'}
                                            )

    custom_forecast_text_final = widgets.BoundedFloatText(
                                                value=0,
                                                min=-1.0,
                                                max=1.0,
                                                step=0.01,
                                                description='Change in county composition:',
                                                disabled=True,
                                                continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,
                                                readout_format='.1f',
                                                style={'description_width': 'initial'}
                                            )

    custom_forecast_text_year = widgets.BoundedFloatText(
                                                value=0,
                                                min=-1.0,
                                                max=1.0,
                                                step=0.01,
                                                description='Change in county composition:',
                                                disabled=True,
                                                continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,
                                                readout_format='.1f',
                                                style={'description_width': 'initial'}
                                            )

    file = open("./red_blue_colorbar.png", "rb");red_blue_colorbar_img = file.read()
    file = open("./purples_colorbar.png", "rb");purples_colorbar_img = file.read()
    colorbar_image = widgets.Image(
        value=purples_colorbar_img,
        format='png',
        width=300,
        height=400,
    )
    race_selector.value = 'white'
    decade_selector.value = '2020'
    year_slider = widgets.IntSlider(min=2010,
                                    max=2100,
                                    step=1,
                                    description='Year',
                                    orientation='horizontal')
    header = widgets.HTML("<h1>Forecasts of neighborhood compositions</h1>", layout=widgets.Layout(height='auto'))
    header.style.text_align = 'center'
    html_out = widgets.HTML(
        value='',
        layout=widgets.Layout(width='auto', height='auto')
    )

    # def get_column_name():
    #     nonlocal decade_selector, race_selector, choropleth_selector
    #     if choropleth_selector.value=='Raw mean compositions': #Choosing just a decade, not the dynamic change
    #         name = "fraction_" + race_selector.value + "_" + decade_selector.value
    #     elif choropleth_selector.value=='Mean change in composition':
    #         name = "delta_fraction_" + race_selector.value + "_" + decade_selector.value[2:4] + \
    #             "_" + decade_selector.value[8:10]
    #     return name

    ####### INITIALIZE MAP ##################################
    geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup, geoid_county)
    m = ipyleaflet.Map(center=initial_position, zoom=10)
    # Choropleth layer for data in column_name
    marker = ipyleaflet.Marker(location=initial_position, draggable=True)
    m.add_control(ipyleaflet.SearchControl(
        position="topleft",
        url='https://nominatim.openstreetmap.org/search?format=json&q={s}',
        zoom=10,
        marker=marker
    ))
    m.add_layer(marker)

    base_layers = m.layers
    column_name = get_column_name()
    layer_blockgroup = ipyleaflet.Choropleth(
        geo_data=geojson_blockgroup,
        choro_data=dict(zip(df_blockgroup_county[:]['geoid_blockgroup'].astype(str), df_blockgroup_county[:][column_name])),
        colormap=linear.Purples_05,
        stroke='False',
        value_min=0,
        value_max=1,
        style={'fillOpacity': 0.8, 'weight': 0})
    m.add_layer(layer_blockgroup)
    m.add_control(decade_selector)

    # Now initialize the line plots
    lin_x = bq.LinearScale()
    lin_y = bq.LinearScale()
    lin_x.min = 0.0
    lin_x.max = 1
    ax_x_bg_forecast = bq.Axis(label='Composition, n', scale=lin_x)
    ax_y_bg_forecast = bq.Axis(label='Probability, P', scale=lin_y, orientation='vertical',label_offset = '70px')
    line_probability = bq.Lines(x=np.linspace(0, 1, len(forecasts[race_selector.value][year_slider.value - year_slider.min, :])),
                           y=forecasts[race_selector.value][year_slider.value - year_slider.min, :],
                           scales={'x': lin_x, 'y': lin_y})
    max_prob = max(forecasts[race_selector.value][year_slider.value - year_slider.min, :])
    n_1990 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_90"]
    n_2000 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_00"]
    n_2010 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_10"]
    line_1990 = bq.Lines(x=[n_1990, n_1990], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted')
    line_2000 = bq.Lines(x=[n_2000, n_2000], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted')
    line_2010 = bq.Lines(x=[n_2010, n_2010], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted')
    lbl_1990 = bq.Label(x=[n_1990], y=[0.25*max_prob], scales={'x': lin_x, 'y': lin_y}, text=['1990'],align = "middle",colors = ['black'])
    lbl_2000 = bq.Label(x=[n_2000], y=[0.5*max_prob], scales={'x': lin_x, 'y': lin_y}, text=['2000'], align="middle",colors = ['black'])
    lbl_2010 = bq.Label(x=[n_2010], y=[0.75*max_prob], scales={'x': lin_x, 'y': lin_y}, text=['2010'], align="middle",colors = ['black'])

    margin_fig = dict(left=100, top=50, bottom=50, right=100)
    fig_bg_forecast = bq.Figure(axes=[ax_x_bg_forecast, ax_y_bg_forecast],
                                marks=[line_probability,line_1990,line_2000,line_2010,lbl_1990,lbl_2000,lbl_2010],
                                fig_margin=margin_fig)

    app = widgets.AppLayout(center=widgets.VBox([decade_selector,m,colorbar_image]),
                            header=widgets.VBox([header,race_selector]),
                            right_sidebar=widgets.VBox([validation_forecast_selector,
                                                        custom_forecast_checkbox,
                                                        custom_forecast_text,
                                                        year_slider,
                                                        fig_bg_forecast]),
                            footer=None,
                            pane_widths=[0, 1, 1],
                            pane_heights=['100px', 4, 0],
                            height='600px',
                            grid_gap="30px")

    #TODO: Add in a interaction for clicking on the map?
    ################## INTERACTIONS ##############################
    # def handle_interaction(**kwargs):
    #     if kwargs['type'] == 'click':
    #         generate_temp_series(*kwargs['coordinates'])
    #         msg = '%s Selected coordinates: %s, Temp: %d C Precipitation: %d mm\n' % (
    #             kwargs['coordinates'], random.randint(-20, 20), random.randint(0, 100))
    #         out.value = add_log(msg)
    #
    # m.on_interaction(handle_interaction)

    def on_location_changed(event):
        location = event['new']
        geoid_county_new = geo_df_county.loc[geo_df_county.contains(Point(location[1], location[0]))]
        geoid_county_new = geoid_county_new.iloc[0]['GEOID']
        nonlocal geo_df_blockgroup, geojson_blockgroup, df_blockgroup, base_layers
        nonlocal line_probability
        nonlocal year_slider, forecasts, geoid_county
        nonlocal geoid_blockgroup, data_selected_blockgroup,n_1990,n_2000,n_2010
        nonlocal initial_position, column_name, years
        nonlocal df_blockgroup, df_blockgroup_county, layer_blockgroup, colorbar_image
        # nonlocal fig_county_forecast, county_heatmap_plots
        if geoid_county != geoid_county_new:
            geoid_county = geoid_county_new
            geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup, geoid_county)
            column_name = get_column_name()
            # county_heatmap_plots = generate_heatmap_plots(df_blockgroup_county, forecasting_parameters, years)
            # fig_county_forecast = county_heatmap_plots[race_selector.value][year_slider.value]
            if column_name[0:5] == "delta":
                value_min = -0.2
                value_max = 0.2
                cmap = linear.RdBu_09
                colorbar_image.value = red_blue_colorbar_img
            else:
                value_min = 0
                value_max = 1
                cmap = linear.Purples_05
                colorbar_image.value = purples_colorbar_img
            layer_blockgroup = ipyleaflet.Choropleth(
                geo_data=geojson_blockgroup,
                choro_data=dict(zip(df_blockgroup_county[:]['geoid_blockgroup'].astype(str), df_blockgroup_county[:][column_name])),
                colormap=cmap,
                stroke='False',
                value_min=value_min,
                value_max=value_max,
                style={'fillOpacity': 0.8, 'weight': 0})
            # all_layers = all_layers + (layer_blockgroup,)
            m.layers = base_layers
            m.add_layer(layer_blockgroup)

            # TODO: Plot the county wide forecast
        geoid_blockgroup = geo_df_blockgroup.loc[
            geo_df_blockgroup.contains(Point(location[1], location[0]))]
        geoid_blockgroup = geoid_blockgroup.iloc[0]['GEOID10']
        data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]
        n_1990 = data_selected_blockgroup.iloc[0]["fraction_" + race_selector.value + "_90"]
        n_2000 = data_selected_blockgroup.iloc[0]["fraction_" + race_selector.value + "_00"]
        n_2010 = data_selected_blockgroup.iloc[0]["fraction_" + race_selector.value + "_10"]
        mus = mus_df[mus_df["geoid_county"] == geoid_county]
        taus = taus_df[taus_df["geoid_county"] == geoid_county]
        hbars = hbars_df[hbars_df["geoid_county"] == geoid_county]
        forecasts = dfft_fcast.generate_forecasts_selected_blockgroup(mus, taus, hbars, data_selected_blockgroup, years)
        # lines_black.y = forecasts["black"][year_slider.value-year_slider.min]
        # lines_hispanic.y = forecasts["hispanic"][year_slider.value - year_slider.min]
        line_probability.y = forecasts["white"][year_slider.value - year_slider.min]
        line_probability.x = np.linspace(0, 1, len(forecasts[race_selector.value][year_slider.value - year_slider.min, :]))
        max_prob = max(forecasts["white"][year_slider.value - year_slider.min])
        line_1990.x = [n_1990,n_1990];line_2000.x = [n_2000,n_2000];line_2010.x = [n_2010,n_2010]
        line_1990.y = [0, max_prob];line_2000.y = [0, max_prob];line_2010.y = [0, max_prob]
        lbl_1990.x = [n_1990];lbl_2000.x = [n_2000];lbl_2010.x = [n_2010]
        lbl_1990.y = [0.25*max_prob];lbl_2000.y = [0.5*max_prob];lbl_2010.y = [0.75*max_prob]
        # vline_white.x = [1, 1]
        return
        # Generate forecasts of the years 2001 to 2019

    def on_race_decade_change(change):
        if change['type'] == 'change' and change['name'] == 'value':
            nonlocal layer_blockgroup, df_blockgroup_county, geojson_blockgroup
            col = get_column_name()
            if col[0:5] == "delta":
                value_min = -0.2
                value_max = 0.2
                cmap = linear.RdBu_09
                colorbar_image.value = red_blue_colorbar_img
            else:
                value_min = 0
                value_max = 1
                cmap = linear.Purples_05
                colorbar_image.value = purples_colorbar_img
            layer_blockgroup = ipyleaflet.Choropleth(
                geo_data=geojson_blockgroup,
                choro_data=dict(zip(df_blockgroup_county[:]['geoid_blockgroup'].astype(str), df_blockgroup_county[:][col])),
                colormap=cmap,
                stroke='False',
                value_min=value_min,
                value_max=value_max,
                style={'fillOpacity': 0.8, 'weight': 0})
            # all_layers = all_layers + (layer_blockgroup,)
            m.layers = base_layers
            m.add_layer(layer_blockgroup)

    def update_year(change):
        if change['type'] == 'change' and change['name'] == 'value':
            line_probability.y = forecasts[race_selector.value][year_slider.value - year_slider.min]
            line_probability.x = np.linspace(0, 1, len(forecasts[race_selector.value][year_slider.value - year_slider.min, :]))
            max_prob = max(forecasts["white"][year_slider.value - year_slider.min])
            line_1990.x = [n_1990, n_1990];line_2000.x = [n_2000, n_2000];line_2010.x = [n_2010, n_2010]
            line_1990.y = [0, max_prob];line_2000.y = [0, max_prob];line_2010.y = [0, max_prob]
            lbl_1990.y = [0.25 * max_prob];lbl_2000.y = [0.5 * max_prob];lbl_2010.y = [0.75 * max_prob]

    # def update_plot_type(change):
    #     print(change)
    #     # if change['type'] == 'change' and change['name'] == 'value':
    #     nonlocal app
    #     if change['new'] == "Neighborhood":
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #     elif change['new'] == "County":
    #         right_sidebar = widgets.VBox([fig_county_forecast, year_slider, plot_type_toggle])
    #     else:
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #
    #     app = widgets.AppLayout(center=m,
    #                             header=header,
    #                             left_sidebar=widgets.VBox([widgets.Label("Race/Ethnicity:"),
    #                                                        race_selector,
    #                                                        widgets.Label("Year:"),
    #                                                        decade_selector]),
    #                             right_sidebar=right_sidebar,
    #                             footer=out,
    #                             pane_widths=['80px', 1, 1],
    #                             pane_heights=['80px', 4, 1],
    #                             height='600px',
    #                             grid_gap="30px")
    #

    ##### INTERACTION FUNCTION CALLS ##############
    #TODO: Add in call to on_location_changed when use search tool for new place

    #TODO: Add in to "on_race_decade_change" code to also update the lines
    race_selector.observe(on_race_decade_change)
    decade_selector.observe(on_race_decade_change)
    year_slider.observe(update_year)
    marker.observe(on_location_changed, 'location')
    # plot_type_toggle.observe(update_plot_type)

    return app


def build_validation_app():
    """
    BRAINSTORMING:
    A user can slider the year slider, change the layer, change the race, or change the marker position.
    1.Changing the year slider only affects the line_vars. Shouldn't need to generate a new forecast or change choropleth
    2.Changing the choropleth layer changes the choropleth but not the line
    3.Changing the race changes the choropleth and the line
    4.Changing the location could change the choropleth and does change the line

    PLANS: A user need only run this lines of code and the user will be presented with
    an interactive map demonstrating both how compositions of cities changed from 2000-2010
    and how our forecasts matched with those observations for White/non-White, Black/non-Black
    and Hispanic/non-Hispanic populations.
    User workflow:
    First, user will need to use the search tool to place the marker within a county of interest.
    Second, user will be presented with an animated heatmap for the probability of observing a
    future composition given the initial composition for a neighborhood of 1000 persons
    Third, user can choose from a dropdown menu to see 9 separate entries for the 2000,2010 observed
    and 2010 forecasted mean for the selected county.
    Fourth, user will be able to click on a neighborhood of interest and see an animated plot for
    the dynamics of that specific neighborhood
    To explore a new county, the user need only use the search tool to find a new county.
    BACKEND:
    First, I need a function that is called when the search tool is used. A county should then be
    selected and the shapefile, Headache, optimum delta_V, and optimum delta_T should be loaded
    from file or calculated. This will be a dataframe.
    Second,

    App widgets...
    1) race_selector: lets user choose whether they want Black/non-Black,White/non-White,Hispanic/non-Hispanic
    2) decade_selector: lets user choose whteher to see just a single decade or the change between two decades
    3) SearchControl: Lets users find new location
    4) marker: Used to collect the location of a given county for analysis
    5) year_slider: Lets users forecast further into the future and intervening years

    :return:
    """
    location = (42.33847, -83.14681) #Somewhere in Wayne county
    years = np.linspace(2000, 2010, 11, 'Int')

    ###### LOADING DATA ################
    geo_df_county = dfft_io.load_county_GIS()
    hbars_df = dfft_io.load_dataframe_with_hbar_from_csv("./parameters/hbar_df.csv")
    hbars_df = hbars_df[hbars_df["decade"]==2000] # hbars for 2000
    mus_df = pd.read_csv("./parameters/mus_df_validation.csv")  # mus used to forecast 2000-2010
    taus_df = pd.read_csv("./parameters/taus_df_90_00.csv")  # taus obtained from 1990-2000 dynamics
    df_blockgroup = pd.read_pickle('./data/forecast_2000_2010_optimized_taus.pickle')

    ######## WIDGETS ############
    # race_selector = widgets.Dropdown(options=('black', 'hispanic', 'white'),
    race_selector = widgets.ToggleButtons(options=(['black','hispanic','white']),
                                          # description = 'Race/Ethnicity:',
                                          # layout=widgets.Layout(width='20%'),
                                          style={'description_width': 'initial'},
                                          disabled=False,
                                          )
    decade_selector = widgets.Dropdown(options=('2000', '2010', '2000->2010 Observed','2000->2010 Forecasted','2000->2010 Error'),
                                       # description='Layer:',
                                       # layout=widgets.Layout(width='20% '),
                                       style={'description_width': 'initial'})
    file = open("./misc/red_blue_colorbar_comp_change.png", "rb");red_blue_colorbar_img_comp_change = file.read()
    file = open("./misc/red_blue_colorbar_DFFT_error.png", "rb");red_blue_colorbar_img_DFFT_error = file.read()
    file = open("./misc/purples_colorbar_composition.png", "rb");purples_colorbar_img_composition = file.read()
    colorbar_image = widgets.Image(
        value=purples_colorbar_img_composition,
        format='png',
        width=300,
        height=40,
    )
    race_selector.value = 'white'
    decade_selector.value = '2000'
    year_slider = widgets.IntSlider(min=min(years),
                                    max=max(years),
                                    step=1,
                                    description='Year',
                                    orientation='horizontal')

    title_text = widgets.HTML(
        value='<h2>Validating DFFT forecasts of 2010 compositions<h2>',
        # layout=widgets.Layout(width='auto', height='auto')
    )

    # header = widgets.HTML("<h1>Forecasts of neighborhood compositions</h1>", layout=widgets.Layout(height='auto'))
    # header.style.text_align = 'center'
    html_out = widgets.HTML(
        value='Loading forecast of initial neighborhood',
        layout=widgets.Layout(width='auto', height='auto')
    )

    def update_text():
        nonlocal html_out
        race = race_selector.value
        colnames_data = ["geoid_blockgroup",race+"_00","total_00",race+"_10","total_10"]
        colnames_data_lbl = ["GEOID 2010", race + " in 2000", "total in 2000", race + " in 2010", "total in 2010"]
        blockgroup_string_data = 'Neighborhood data:'
        for colname,lbl in zip(colnames_data,colnames_data_lbl):
            blockgroup_string_data += ' ' + lbl + ': ' + str(round(data_selected_blockgroup[colname].iloc[0],3)) + ','
        colnames_forecast = [race+"_2010_mean",race+"_2010_variance",race+"_2010_skewness"]
        colnames_forecast_lbl = [race + " forecasted mean 2010", race + " forecasted variance 2010",race + " forecasted skewness 2010"]
        blockgroup_string_forecast = 'Forecast results:'
        for colname,lbl in zip(colnames_forecast,colnames_forecast_lbl):
            blockgroup_string_forecast += ' ' + lbl + ': ' + str(round(data_selected_blockgroup[colname].iloc[0],3)) + ','
        html_out.value = blockgroup_string_data + '<br />' + blockgroup_string_forecast
        return

    def thinking_text():
        nonlocal html_out
        html_out.value = "THINKING"
        return

    def get_column_name():
        nonlocal decade_selector, race_selector
        race = race_selector.value
        layer_type = decade_selector.value
        if layer_type=='2000':
            name = "fraction_" + race + "_" + layer_type
        elif layer_type == '2010':  # Choosing just a decade, not the dynamic change
            name = "fraction_" + race + "_" + layer_type
        elif layer_type == '2000->2010 Observed':  # Choosing just a decade, not the dynamic change
            name = "delta_fraction_observed_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Forecasted':  #Choosing just a decade, not the dynamic change
            name = "delta_fraction_forecasted_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Error':  # Choosing just a decade, not the dynamic change
            name = "error_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Log-likelihood':  # Choosing just a decade, not the dynamic change
            name = "lglik_" + race + "_2000_2010"
        return name

    ####### INITIALIZE MAP ##################################
    geoid_county = geo_df_county.loc[geo_df_county.contains(Point(location[1], location[0]))]
    geoid_county = int(geoid_county.iloc[0]['GEOID'])
    # print(mus_df)
    # mus = mus_df[mus_df["geoid_county"] == geoid_county]
    # print(mus)
    # taus = taus_df[taus_df["geoid_county"] == geoid_county]
    # hbars = hbars_df[hbars_df["county"] == geoid_county]
    geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup, geoid_county)
    m = ipyleaflet.Map(center=location, zoom=10)
    # Choropleth layer for data in column_name
    marker = ipyleaflet.Marker(location=location, draggable=True)
    m.add_control(ipyleaflet.SearchControl(
        position="topleft",
        url='https://nominatim.openstreetmap.org/search?format=json&q={s}',
        zoom=10,
        marker=marker
    ))
    m.add_layer(marker)

    base_layers = m.layers
    column_name = get_column_name()

    def update_choropleth():
        nonlocal column_name, geojson_blockgroup,df_blockgroup_county,m,base_layers
        column_name = get_column_name()
        if column_name[0:5] == "delta":
            value_min = -0.1; value_max = 0.1; cmap = linear.RdBu_09
            colorbar_image.value = red_blue_colorbar_img_comp_change
        elif column_name[0:5] == "error":
            value_min = -3; value_max = 3; cmap = linear.RdBu_09
            colorbar_image.value = red_blue_colorbar_img_DFFT_error
        elif column_name[0:8] == "fraction":
            value_min = 0;value_max = 1;cmap = linear.Purples_05
            colorbar_image.value = purples_colorbar_img_composition
        layer_blockgroup = ipyleaflet.Choropleth(
            geo_data=geojson_blockgroup,
            choro_data=dict(zip(df_blockgroup_county[:]['geoid_blockgroup'].astype(str), df_blockgroup_county[:][column_name])),
            colormap=cmap,
            stroke='False',
            value_min=value_min,
            value_max=value_max,
            style={'fillOpacity': 0.8, 'weight': 0})
        # all_layers = all_layers + (layer_blockgroup,)
        m.layers = base_layers
        m.add_layer(layer_blockgroup)
        return
    update_choropleth()

    ########## INITIAL NEIGHBORHOOD LINE PLOT ###################
    # Now initialize the line plots and collect into a dictionary
    lin_x = bq.LinearScale()
    lin_y = bq.LinearScale()
    lin_x.min = 0.0
    lin_x.max = 1
    ax_x_bg_forecast = bq.Axis(label='Composition, n', scale=lin_x,tick_values = np.array([0.0,0.2,0.4,0.6,0.8,1.0]))
    ax_y_bg_forecast = bq.Axis(label='Probability, P', scale=lin_y, orientation='vertical',label_offset = '70px')
    line_vars = (bq.Lines(x=np.linspace(0, 1, 1001),y=np.zeros(1001),scales={'x': lin_x, 'y': lin_y}),
                 bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 bq.Label(x=[0], y=[0.25], scales={'x': lin_x, 'y': lin_y}, text=['1990'],align = "middle",colors = ['black']),
                 bq.Label(x=[0], y=[0.5], scales={'x': lin_x, 'y': lin_y}, text=['2000'], align="middle",colors = ['black']),
                 bq.Label(x=[0], y=[0.75], scales={'x': lin_x, 'y': lin_y}, text=['2010'], align="middle",colors = ['black'])
                )

    geoid_blockgroup = geo_df_blockgroup.loc[
        geo_df_blockgroup.contains(Point(location[1], location[0]))]
    geoid_blockgroup = int(geoid_blockgroup.iloc[0]['GEOID10'])
    data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]
    forecasts = dfft_fcast.generate_blockgroup_forecast(race_selector.value, mus_df[mus_df["geoid_county"] == geoid_county],
                                                        taus_df[taus_df["geoid_county"] == geoid_county],
                                                        hbars_df[hbars_df["county"] == geoid_county],
                                                        data_selected_blockgroup, years)

    def update_forecast_lines():
        nonlocal location
        # Should be called whenever the race is changed or the location is changed.
        nonlocal geoid_county,line_vars,race_selector
        nonlocal geo_df_blockgroup,geojson_blockgroup,data_selected_blockgroup, forecasts, geoid_blockgroup
        line_probability,line_1990,line_2000,line_2010,lbl_1990,lbl_2000,lbl_2010 = line_vars
        # Find the blockgroup that matches the location
        # Generate forecasts of that neighborhood
        geoid_blockgroup = geo_df_blockgroup.loc[geo_df_blockgroup.contains(Point(location[1], location[0]))]
        geoid_blockgroup = int(geoid_blockgroup.iloc[0]['GEOID10'])
        data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]
        forecasts = dfft_fcast.generate_blockgroup_forecast(race_selector.value,mus_df[mus_df["geoid_county"]==geoid_county],
                                                            taus_df[taus_df["geoid_county"]==geoid_county],
                                                            hbars_df[hbars_df["county"]==geoid_county],
                                                            data_selected_blockgroup, years)
        max_prob = max(forecasts[year_slider.value - year_slider.min, :])
        n_1990 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_1990"]
        n_2000 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_2000"]
        n_2010 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_2010"]
        line_probability.y = forecasts[year_slider.value - year_slider.min]
        np.save("tmp",forecasts[9])
        np.save("tmp_2000", n_2000)
        line_probability.x = np.linspace(0, 1, len(forecasts[year_slider.value - year_slider.min,:]))
        line_1990.x = [n_1990,n_1990];line_2000.x = [n_2000,n_2000];line_2010.x = [n_2010,n_2010]
        line_1990.y = [0, max_prob];line_2000.y = [0, max_prob];line_2010.y = [0, max_prob]
        lbl_1990.x = [n_1990];lbl_2000.x = [n_2000];lbl_2010.x = [n_2010]
        lbl_1990.y = [0.25*max_prob];lbl_2000.y = [0.5*max_prob];lbl_2010.y = [0.75*max_prob]
        return forecasts

    update_forecast_lines()
    update_text()
    line_probability, line_1990, line_2000, line_2010, lbl_1990, lbl_2000, lbl_2010 = line_vars

    ########## BUILD APP ##########################
    margin_fig = dict(left=100, top=50, bottom=50, right=100)
    fig_bg_forecast = bq.Figure(axes=[ax_x_bg_forecast, ax_y_bg_forecast],
                                marks=[line_vars[0],line_vars[1],line_vars[2],line_vars[3],line_vars[4],line_vars[5],line_vars[6]],
                                fig_margin=margin_fig)


    tmpbox = widgets.HBox([race_selector,decade_selector])

    app = widgets.AppLayout(
                            left_sidebar=None,
                            center=widgets.VBox([m,colorbar_image]),
                            header=widgets.VBox([title_text, tmpbox]),
                            # header=widgets.HBox([race_selector,decade_selector]),
                            # header=race_selector,
                            right_sidebar=widgets.VBox([fig_bg_forecast,year_slider]),
                            footer=html_out,
                            pane_widths=['0px','400px', '370px'],
                            # pane_widths=['160px', '400px', '370px'],
                            pane_heights=['120px', '300px', '100px'],
                            height='520px',
                            width='800px',
                            # grid_gap="30px"
                            )

    # app = widgets.AppLayout(left_sidebar=widgets.VBox([race_selector,decade_selector]),
    #                         center=m,
    #                         # header=widgets.VBox([race_selector]),
    #                         right_sidebar=widgets.VBox([fig_bg_forecast,year_slider]),
    #                         footer=html_out,
    #                         pane_widths=['160px', '400px', '370px'],
    #                         # pane_heights=['100px', 4, 0],
    #                         height='370px',
    #                         width='1000px',
    #                         # grid_gap="30px"
    #                         )

    ################## INTERACTIONS ##############################
    def on_location_changed(event):
        # Steps
        # 1. Check to see if the county is changed. If so, update the choropleth and county specific variables
        # 2. Generate forecasts update the line plots
        thinking_text()
        nonlocal location
        location = event['new']
        nonlocal geoid_county, html_out, geojson_blockgroup,df_blockgroup_county,geo_df_blockgroup,df_blockgroup
        # nonlocal mus,taus,hbars
        geoid_county_new = geo_df_county.loc[geo_df_county.contains(Point(location[1], location[0]))]
        if geoid_county_new.size==0:
            html_out.value = 'Marker not found within a US county'
            return
        geoid_county_new = int(geoid_county_new.iloc[0]['GEOID'])
        if geoid_county != geoid_county_new:
            geoid_county = geoid_county_new
            geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup,
                                                                                             geoid_county)
            update_choropleth()
        update_forecast_lines()
        update_text()
        return

    def on_race_decade_change(change):
        #1. Update the choropleth
        #2. Update the forecast lines
        if change['type'] == 'change' and change['name'] == 'value':
            thinking_text()
            update_choropleth()
            update_forecast_lines()
            update_text()
        return

    def update_year(change):
        # nonlocal race_selector
        # race = race_selector.value
        line_probability, line_1990, line_2000, line_2010, lbl_1990, lbl_2000, lbl_2010 = line_vars
        if change['type'] == 'change' and change['name'] == 'value':
            line_probability.y = forecasts[year_slider.value - year_slider.min]
            line_probability.x = np.linspace(0, 1, len(forecasts[year_slider.value - year_slider.min, :]))
            max_prob = max(forecasts[year_slider.value - year_slider.min])
            line_1990.y = [0, max_prob];line_2000.y = [0, max_prob];line_2010.y = [0, max_prob]
            lbl_1990.y = [0.25 * max_prob];lbl_2000.y = [0.5 * max_prob];lbl_2010.y = [0.75 * max_prob]

    # def update_plot_type(change):
    #     print(change)
    #     # if change['type'] == 'change' and change['name'] == 'value':
    #     nonlocal app
    #     if change['new'] == "Neighborhood":
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #     elif change['new'] == "County":
    #         right_sidebar = widgets.VBox([fig_county_forecast, year_slider, plot_type_toggle])
    #     else:
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #
    #     app = widgets.AppLayout(center=m,
    #                             header=header,
    #                             left_sidebar=widgets.VBox([widgets.Label("Race/Ethnicity:"),
    #                                                        race_selector,
    #                                                        widgets.Label("Year:"),
    #                                                        decade_selector]),
    #                             right_sidebar=right_sidebar,
    #                             footer=out,
    #                             pane_widths=['80px', 1, 1],
    #                             pane_heights=['80px', 4, 1],
    #                             height='600px',
    #                             grid_gap="30px")
    #

    ##### INTERACTION FUNCTION CALLS ##############
    #TODO: Add in call to on_location_changed when use search tool for new place

    #TODO: Add in to "on_race_decade_change" code to also update the lines
    race_selector.observe(on_race_decade_change)
    decade_selector.observe(on_race_decade_change)
    year_slider.observe(update_year)
    marker.observe(on_location_changed, 'location')
    # plot_type_toggle.observe(update_plot_type)

    return app


def build_forecasting_app():
    """
    BRAINSTORMING:
    A user can slider the year slider, change the layer, change the race, or change the marker position.
    1.Changing the year slider only affects the line_vars. Shouldn't need to generate a new forecast or change choropleth
    2.Changing the choropleth layer changes the choropleth but not the line
    3.Changing the race changes the choropleth and the line
    4.Changing the location could change the choropleth and does change the line

    PLANS: A user need only run this lines of code and the user will be presented with
    an interactive map demonstrating both how compositions of cities changed from 2000-2010
    and how our forecasts matched with those observations for White/non-White, Black/non-Black
    and Hispanic/non-Hispanic populations.
    User workflow:
    First, user will need to use the search tool to place the marker within a county of interest.
    Second, user will be presented with an animated heatmap for the probability of observing a
    future composition given the initial composition for a neighborhood of 1000 persons
    Third, user can choose from a dropdown menu to see 9 separate entries for the 2000,2010 observed
    and 2010 forecasted mean for the selected county.
    Fourth, user will be able to click on a neighborhood of interest and see an animated plot for
    the dynamics of that specific neighborhood
    To explore a new county, the user need only use the search tool to find a new county.
    BACKEND:
    First, I need a function that is called when the search tool is used. A county should then be
    selected and the shapefile, Headache, optimum delta_V, and optimum delta_T should be loaded
    from file or calculated. This will be a dataframe.
    Second,

    App widgets...
    1) race_selector: lets user choose whether they want Black/non-Black,White/non-White,Hispanic/non-Hispanic
    2) decade_selector: lets user choose whteher to see just a single decade or the change between two decades
    3) SearchControl: Lets users find new location
    4) marker: Used to collect the location of a given county for analysis
    5) year_slider: Lets users forecast further into the future and intervening years

    :return:
    """
    location = (42.33847, -83.14681) #Somewhere in Wayne county
    years = np.linspace(2010, 2050, 9, 'Int')

    ###### LOADING DATA ################
    geo_df_county = dfft_io.load_county_GIS()
    # TODO: WARNING, Add in the 2010 hbars to hbar.csv, right now using the wrong hbars for forecasting app
    hbars_df = dfft_io.load_dataframe_with_hbar_from_csv("./parameters/hbar_df.csv")
    hbars_df = hbars_df[hbars_df["decade"]==2000] # hbars for 2000
    # TODO: WARNING, need to use correct taus and eventually load in df_blockgroup with SSP forecasts. For now, I just do the custom forecast

    # mus_df = pd.read_csv("./parameters/mus_df_SSP2.csv")  # mus used to forecast 2020-2100. Also used for the forecasts
    # taus_df = pd.read_csv("./parameters/taus_df_forecasting.csv")  # taus obtained from 2000-2010 dynamics
    # df_blockgroup = pd.read_csv("./data/forecasted_moments_SSP2_df.csv")  # each block group mean 2020-2100
    mus_df = pd.read_csv("./parameters/mus_df_validation.csv")  # mus used to forecast 2000-2010
    taus_df = pd.read_csv("./parameters/taus_df_00_10.csv")  # taus obtained from 1990-2000 dynamics
    df_blockgroup = pd.read_pickle('./data/forecast_2000_2010_optimized_taus.pickle')

    ######## WIDGETS ############
    race_selector = widgets.ToggleButtons(options=(['black','hispanic','white']),
                                          # description = 'Race/Ethnicity:',
                                          style={'description_width': 'initial'},
                                          disabled=False,)
    decade_selector = widgets.Dropdown(options=('2010',),#, '2010->2020','2010->2020','2010->2030','2010->2040','2010->2050','2010->2060','2010->2070','2010->2080','2010->2090','2010->2100'),
                                       # description='Layer:',
                                       # layout=widgets.Layout(width='20% '),
                                       style={'description_width': 'initial'})
    projection_selector = widgets.Dropdown(options=('SSP1', 'SSP2', 'SSP3', 'SSP4', 'SSP5'),
                                           description='County projection:',
                                           # layout=widgets.Layout(width='20% '),
                                           style={'description_width': 'initial'},
                                           disabled=True,
                                           layout=widgets.Layout(width='200px'))
    custom_forecast_checkbox = widgets.Checkbox(value=True,
                                        description='Custom county forecast',
                                        disabled=True,
                                        indent=False,
                                        layout=widgets.Layout(width='200px'))
    custom_forecast_text_initial = widgets.BoundedFloatText(value=0.0,min=0.0,max=1.0,
                                                step=0.01,description='2010 Composition:',
                                                disabled=True,continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,readout_format='.1f',
                                                style={'description_width': 'initial'},
                                                layout=widgets.Layout(width='200px'))

    custom_forecast_text_final = widgets.BoundedFloatText(value=0,min=0.0,max=1.0,step=0.01,
                                                description='2050 Composition:',
                                                disabled=False,
                                                continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,
                                                readout_format='.1f',
                                                style={'description_width': 'initial'},
                                                layout=widgets.Layout(width='200px'))

    custom_forecast_text_year = widgets.BoundedFloatText(
                                                value=0,
                                                min=-1.0,
                                                max=1.0,
                                                step=0.01,
                                                description='Change in county composition:',
                                                disabled=True,
                                                continuous_update=False,
                                                orientation='horizontal',
                                                readout=True,
                                                readout_format='.1f',
                                                style={'description_width': 'initial'}
                                            )
    file = open("./misc/red_blue_colorbar_comp_change.png", "rb");red_blue_colorbar_img_comp_change = file.read()
    file = open("./misc/red_blue_colorbar_DFFT_error.png", "rb");red_blue_colorbar_img_DFFT_error = file.read()
    file = open("./misc/purples_colorbar_composition.png", "rb");purples_colorbar_img_composition = file.read()
    colorbar_image = widgets.Image(
        value=purples_colorbar_img_composition,
        format='png',
        width=300,
        height=40,
    )
    race_selector.value = 'white'
    decade_selector.value = '2010'
    year_slider = widgets.IntSlider(min=min(years),max=max(years),
                                    step=years[1]-years[0],description='Year',orientation='horizontal')
    header = widgets.HTML("<h1>Forecasts of neighborhood compositions</h1>", layout=widgets.Layout(height='auto'))
    header.style.text_align = 'center'
    html_out = widgets.HTML(
        value='Loading forecast of initial neighborhood',
        layout=widgets.Layout(width='auto', height='auto')
    )
    title_text = widgets.HTML(value='<h2>Forecasting neighborhood compositions into the future<h2>',)
    projection_title_text = widgets.HTML(value='<h4>County scale projection<h4>',)

    def update_text():
        nonlocal html_out
        race = race_selector.value
        colnames_data = ["geoid_blockgroup",race+"_10","total_10"]
        colnames_data_lbl = ["GEOID 2010", race + " in 2010", "total in 2010"]
        blockgroup_string_data = 'Neighborhood data:'
        for colname,lbl in zip(colnames_data,colnames_data_lbl):
            blockgroup_string_data += ' ' + lbl + ': ' + str(round(data_selected_blockgroup[colname].iloc[0],3)) + ','
        # colnames_forecast = [race+"_2010_mean",race+"_2010_variance",race+"_2010_skewness"]
        # colnames_forecast_lbl = [race + " forecasted mean 2010", race + " forecasted variance 2010",race + " forecasted skewness 2010"]
        # blockgroup_string_forecast = 'Forecast results:'
        # for colname,lbl in zip(colnames_forecast,colnames_forecast_lbl):
        #     blockgroup_string_forecast += ' ' + lbl + ': ' + str(round(data_selected_blockgroup[colname].iloc[0],3)) + ','
        # html_out.value = blockgroup_string_data + '<br />' + blockgroup_string_forecast
        html_out.value = blockgroup_string_data

        return

    def thinking_text():
        nonlocal html_out
        html_out.value = "THINKING"
        return

    def get_column_name():
        nonlocal decade_selector, race_selector
        race = race_selector.value
        layer_type = decade_selector.value
        if layer_type=='2000':
            name = "fraction_" + race + "_" + layer_type
        elif layer_type == '2010':  # Choosing just a decade, not the dynamic change
            name = "fraction_" + race + "_" + layer_type
        elif layer_type == '2000->2010 Observed':  # Choosing just a decade, not the dynamic change
            name = "delta_fraction_observed_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Forecasted':  #Choosing just a decade, not the dynamic change
            name = "delta_fraction_forecasted_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Error':  # Choosing just a decade, not the dynamic change
            name = "error_" + race + "_2000_2010"
        elif layer_type == '2000->2010 Log-likelihood':  # Choosing just a decade, not the dynamic change
            name = "lglik_" + race + "_2000_2010"
        return name

    ####### INITIALIZE MAP ##################################
    geoid_county = geo_df_county.loc[geo_df_county.contains(Point(location[1], location[0]))]
    geoid_county = int(geoid_county.iloc[0]['GEOID'])
    # print(mus_df)
    # mus = mus_df[mus_df["geoid_county"] == geoid_county]
    # print(mus)
    # taus = taus_df[taus_df["geoid_county"] == geoid_county]
    # hbars = hbars_df[hbars_df["county"] == geoid_county]
    geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup, geoid_county)
    m = ipyleaflet.Map(center=location, zoom=10)
    # Choropleth layer for data in column_name
    marker = ipyleaflet.Marker(location=location, draggable=True)
    m.add_control(ipyleaflet.SearchControl(
        position="topleft",
        url='https://nominatim.openstreetmap.org/search?format=json&q={s}',
        zoom=10,
        marker=marker
    ))
    m.add_layer(marker)

    base_layers = m.layers
    column_name = get_column_name()
    mu_custom = 0.0

    def update_mu_custom():
        nonlocal mu_custom
        # If the custom checkbox is selected, then when the county changes or the custom forecast value is change
        # optimize a new value for mu.
        taus = taus_df[taus_df["geoid_county"] == geoid_county]
        tau = taus[taus["race"] == race_selector.value]["tau"].iloc[0]
        mu_custom = dfft_fcast.estimate_mu(race_selector.value,df_blockgroup_county,tau,years[-1]-years[0],
                                           hbars_df[hbars_df["county"] == geoid_county],
                                           custom_forecast_text_final.value,
                                           mu_lims = [-1, 1])
        return

    def update_choropleth():
        nonlocal column_name, geojson_blockgroup,df_blockgroup_county,m
        n_county = round(df_blockgroup_county[race_selector.value + '_10'].sum() / df_blockgroup_county['total_10'].sum(),4)
        custom_forecast_text_initial.value = n_county
        custom_forecast_text_final.value = n_county
        column_name = get_column_name()
        if column_name[0:5] == "delta":
            value_min = -0.1; value_max = 0.1; cmap = linear.RdBu_09
            colorbar_image.value = red_blue_colorbar_img_comp_change
        elif column_name[0:5] == "error":
            value_min = -3; value_max = 3; cmap = linear.RdBu_09
            colorbar_image.value = red_blue_colorbar_img_DFFT_error
        elif column_name[0:8] == "fraction":
            value_min = 0;value_max = 1;cmap = linear.Purples_05
            colorbar_image.value = purples_colorbar_img_composition
        layer_blockgroup = ipyleaflet.Choropleth(
            geo_data=geojson_blockgroup,
            choro_data=dict(zip(df_blockgroup_county[:]['geoid_blockgroup'].astype(str), df_blockgroup_county[:][column_name])),
            colormap=cmap,
            stroke='False',
            value_min=value_min,
            value_max=value_max,
            style={'fillOpacity': 0.8, 'weight': 0})
        # all_layers = all_layers + (layer_blockgroup,)
        m.layers = base_layers
        m.add_layer(layer_blockgroup)
        return
    update_choropleth()

    ########## INITIAL NEIGHBORHOOD LINE PLOT ###################
    # Now initialize the line plots and collect into a dictionary
    lin_x = bq.LinearScale()
    lin_y = bq.LinearScale()
    lin_x.min = 0.0
    lin_x.max = 1
    ax_x_bg_forecast = bq.Axis(label='Composition, n', scale=lin_x,tick_values = np.array([0.0,0.2,0.4,0.6,0.8,1.0]))
    ax_y_bg_forecast = bq.Axis(label='Probability, P', scale=lin_y, orientation='vertical',label_offset = '70px')
    line_vars = (bq.Lines(x=np.linspace(0, 1, 1001),y=np.zeros(1001),scales={'x': lin_x, 'y': lin_y}),
                 # bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 # bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 bq.Lines(x=[0, 0], y=[0, 1], scales={'x': lin_x, 'y': lin_y},preserve_domain={'x': True, 'y': False},colors = ['black'], line_style = 'dotted'),
                 # bq.Label(x=[0], y=[0.25], scales={'x': lin_x, 'y': lin_y}, text=['1990'],align = "middle",colors = ['black']),
                 # bq.Label(x=[0], y=[0.5], scales={'x': lin_x, 'y': lin_y}, text=['2000'], align="middle",colors = ['black']),
                 bq.Label(x=[0], y=[0.75], scales={'x': lin_x, 'y': lin_y}, text=['2010'], align="middle",colors = ['black'])
                )

    geoid_blockgroup = geo_df_blockgroup.loc[
        geo_df_blockgroup.contains(Point(location[1], location[0]))]
    geoid_blockgroup = int(geoid_blockgroup.iloc[0]['GEOID10'])
    data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]
    forecasts = dfft_fcast.generate_blockgroup_forecast(race_selector.value, mus_df[mus_df["geoid_county"] == geoid_county],
                                                        taus_df[taus_df["geoid_county"] == geoid_county],
                                                        hbars_df[hbars_df["county"] == geoid_county],
                                                        data_selected_blockgroup, years,
                                                        initial_year=2010)

    def update_forecast_lines():
        nonlocal location
        # Should be called whenever the race is changed or the location is changed.
        nonlocal geoid_county,line_vars,race_selector
        nonlocal geo_df_blockgroup,geojson_blockgroup,data_selected_blockgroup, forecasts, geoid_blockgroup
        # line_probability,line_1990,line_2000,line_2010,lbl_1990,lbl_2000,lbl_2010 = line_vars
        line_probability, line_2010, lbl_2010 = line_vars
        # Find the blockgroup that matches the location
        # Generate forecasts of that neighborhood
        geoid_blockgroup = geo_df_blockgroup.loc[geo_df_blockgroup.contains(Point(location[1], location[0]))]
        geoid_blockgroup = int(geoid_blockgroup.iloc[0]['GEOID10'])
        data_selected_blockgroup = df_blockgroup[df_blockgroup['geoid_blockgroup'] == geoid_blockgroup]
        forecasts = dfft_fcast.generate_blockgroup_forecast(race_selector.value,mu_custom,
                                                            taus_df[taus_df["geoid_county"]==geoid_county],
                                                            hbars_df[hbars_df["county"]==geoid_county],
                                                            data_selected_blockgroup, years,
                                                            initial_year=2010)
        max_prob = max(forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0])), :])
        n_2010 = data_selected_blockgroup.iloc[0]["fraction_"+race_selector.value+"_2010"]
        line_probability.y = forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0]))]
        np.save("tmp",forecasts[8])
        np.save("tmp_2010", n_2010)
        line_probability.x = np.linspace(0, 1, len(forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0])),:]))
        line_2010.x = [n_2010,n_2010]; line_2010.y = [0, max_prob]
        lbl_2010.x = [n_2010]; lbl_2010.y = [0.75*max_prob]
        return forecasts

    update_forecast_lines()
    update_text()

    ########## BUILD APP ##########################
    margin_fig = dict(left=100, top=50, bottom=50, right=100)
    fig_bg_forecast = bq.Figure(axes=[ax_x_bg_forecast, ax_y_bg_forecast],
                                marks=[line_vars[0],line_vars[1],line_vars[2]],
                                fig_margin=margin_fig)

    tmpbox = widgets.HBox([race_selector,decade_selector])
    app = widgets.AppLayout(left_sidebar=widgets.VBox([m,colorbar_image]),
                            # center=widgets.VBox([widgets.HBox([projection_selector,custom_forecast_checkbox]),
                            #                      widgets.HBox([custom_forecast_text_initial, custom_forecast_text_final])]),
                            center=widgets.VBox([projection_title_text, projection_selector, custom_forecast_checkbox,custom_forecast_text_initial, custom_forecast_text_final]),
                            header=widgets.VBox([title_text,tmpbox]),
                            # header=race_selector,
                            right_sidebar=widgets.VBox([fig_bg_forecast,year_slider]),
                            footer=html_out,
                            pane_widths=['345px','210px', '345px'],
                            # pane_widths=['160px', '400px', '370px'],
                            pane_heights=['120px', '300px', '100px'],
                            height='520px',
                            width='900px',
                            # grid_gap="30px"
                            )

    # app = widgets.AppLayout(
    #                         left_sidebar=None,
    #                         center=widgets.VBox([m,colorbar_image]),
    #                         header=widgets.HBox([race_selector,decade_selector]),
    #                         # header=race_selector,
    #                         right_sidebar=widgets.VBox([fig_bg_forecast,year_slider]),
    #                         footer=html_out,
    #                         pane_widths=['0px','400px', '370px'],
    #                         # pane_widths=['160px', '400px', '370px'],
    #                         pane_heights=['50px', '300px', '100px'],
    #                         height='450px',
    #                         width='800px',
    #                         # grid_gap="30px"
    #                         )

    # app = widgets.AppLayout(left_sidebar=widgets.VBox([race_selector,decade_selector]),
    #                         center=m,
    #                         # header=widgets.VBox([race_selector]),
    #                         right_sidebar=widgets.VBox([fig_bg_forecast,year_slider]),
    #                         footer=html_out,
    #                         pane_widths=['160px', '400px', '370px'],
    #                         # pane_heights=['100px', 4, 0],
    #                         height='370px',
    #                         width='1000px',
    #                         # grid_gap="30px"
    #                         )

    ################## INTERACTIONS ##############################
    def on_location_changed(event):
        # Steps
        # 1. Check to see if the county is changed. If so, update the choropleth and county specific variables
        # 2. Generate forecasts update the line plots
        thinking_text()
        nonlocal location
        location = event['new']
        nonlocal geoid_county, html_out, geojson_blockgroup,df_blockgroup_county,geo_df_blockgroup,df_blockgroup
        # nonlocal mus,taus,hbars
        geoid_county_new = geo_df_county.loc[geo_df_county.contains(Point(location[1], location[0]))]
        if geoid_county_new.size==0:
            html_out.value = 'Marker not found within a US county'
            return
        geoid_county_new = int(geoid_county_new.iloc[0]['GEOID'])
        if geoid_county != geoid_county_new:
            geoid_county = geoid_county_new
            geojson_blockgroup, df_blockgroup_county, geo_df_blockgroup = update_county_data(df_blockgroup,
                                                                                             geoid_county)
            update_choropleth()
        update_forecast_lines()
        update_text()
        return

    def on_race_decade_change(change):
        #1. Update the choropleth
        #2. Update the forecast lines
        if change['type'] == 'change' and change['name'] == 'value':
            thinking_text()
            update_choropleth()
            update_forecast_lines()
            update_text()
        return

    def update_year(change):
        # nonlocal race_selector
        # race = race_selector.value
        line_probability, line_2010, lbl_2010 = line_vars
        if change['type'] == 'change' and change['name'] == 'value':
            line_probability.y = forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0])),:]
            line_probability.x = np.linspace(0, 1, len(forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0])), :]))
            max_prob = max(forecasts[int((year_slider.value - year_slider.min)/(years[1]-years[0])),:])
            line_2010.y = [0, max_prob]
            lbl_2010.y = [0.75 * max_prob]

    def on_custom_changed(change):
        if change['type'] == 'change' and change['name'] == 'value':
            thinking_text()
            update_mu_custom()
            update_forecast_lines()
            update_text()

    # def update_plot_type(change):
    #     print(change)
    #     # if change['type'] == 'change' and change['name'] == 'value':
    #     nonlocal app
    #     if change['new'] == "Neighborhood":
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #     elif change['new'] == "County":
    #         right_sidebar = widgets.VBox([fig_county_forecast, year_slider, plot_type_toggle])
    #     else:
    #         right_sidebar = widgets.VBox([fig_bg_forecast, year_slider, plot_type_toggle])
    #
    #     app = widgets.AppLayout(center=m,
    #                             header=header,
    #                             left_sidebar=widgets.VBox([widgets.Label("Race/Ethnicity:"),
    #                                                        race_selector,
    #                                                        widgets.Label("Year:"),
    #                                                        decade_selector]),
    #                             right_sidebar=right_sidebar,
    #                             footer=out,
    #                             pane_widths=['80px', 1, 1],
    #                             pane_heights=['80px', 4, 1],
    #                             height='600px',
    #                             grid_gap="30px")
    #

    ##### INTERACTION FUNCTION CALLS ##############
    #TODO: Add in call to on_location_changed when use search tool for new place

    #TODO: Add in to "on_race_decade_change" code to also update the lines
    race_selector.observe(on_race_decade_change)
    decade_selector.observe(on_race_decade_change)
    year_slider.observe(update_year)
    marker.observe(on_location_changed, 'location')
    custom_forecast_text_final.observe(on_custom_changed)
    # plot_type_toggle.observe(update_plot_type)

    return app


#
# def generate_forecasts_selected_blockgroup(H_dict, steps_per_year_all_races, mus, data_selected_blockgroup, years):
#     """
#     Generate forecasts for the probability of observing a given fraction of persons for binary White, Black, and Hispanic
#     classifications
#     :param H_functions:
#     :param steps_per_year_all_races:
#     :param mus:
#     :param data_selected_blockgroup:
#     :return:
#     """
#     if years[0] == 2000: # Leaving this in for now to allow flexibility later to produce forecasts from 1990 instead
#         total = data_selected_blockgroup['total_00']
#     # elif years[0] == 1990:
#     #     total = data_selected_blockgroup['total_90']
#     transition_matrices = {
#         'black': dfft_fcast.build_transition_matrix(dfft_fcast.shift_H_by_mu(H_dict['black_00'], mus['black_00_10']), total),
#         'hispanic': dfft_fcast.build_transition_matrix(dfft_fcast.shift_H_by_mu(H_dict['hispanic_00'], mus['hispanic_00_10']), total),
#         'white': dfft_fcast.build_transition_matrix(dfft_fcast.shift_H_by_mu(H_dict['white_00'], mus['white_00_10']), total)}
#
#     def forecast_from_transition_matrix(transition_matrix, N_initial, number_years, tot, steps_per_year):
#         state_vector = np.zeros((tot + 1, 1))
#         state_vector[N_initial] = 1
#         transition_matrix = np.linalg.matrix_power(transition_matrix, round(steps_per_year * tot))
#         forecast = np.zeros(tot + 1, number_years)
#         forecast[:, 0] = state_vector
#         for i in range(number_years - 1):
#             state_vector = np.matmul(transition_matrix, state_vector)
#             forecast[:, i + 1] = state_vector
#         return forecast
#
#     forecasts = {'black': forecast_from_transition_matrix(transition_matrices['black'],
#                                                           round(data_selected_blockgroup['black_00']),
#                                                           len(years), total, steps_per_year_all_races['black']),
#                  'hispanic': forecast_from_transition_matrix(transition_matrices['hispanic'],
#                                                              round(data_selected_blockgroup['hispanic_00']),
#                                                              len(years), total, steps_per_year_all_races['hispanic']),
#                  'white': forecast_from_transition_matrix(transition_matrices['white'],
#                                                           round(data_selected_blockgroup['white_00']),
#                                                           len(years), total, steps_per_year_all_races['white'])
#                  }
#     return forecasts
#
#
# rows = []
# h = np.load('dict_h.npy', allow_pickle = True).item()
#
# # appending rows
# for decade_key, h_decade in h.items():
#     for race_key,h_race in h_decade.items():
#         for county_key,h_county in h_race.items():
#             data_row = data['Student']
#             time = data['Name']
#             row = {"decade": decade_key, "race": race_key, "county": county_key, "H": h_county}
#             rows.append(row)
#             print(county_key)
#
# df = pd.DataFrame(rows)


# [
# "Paired", "Set3", "Pastel1", "Set1",
# "Greys", "Greens", "Reds", "Purples", "Oranges", "Blues",
# "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "RdPu",
# "PuRd", "PuBuGn", "PuBu", "OrRd", "GnBu", "BuPu",
# "BuGn", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy",
# "RdYlBu", "RdYlGn", "Spectral"
# ]