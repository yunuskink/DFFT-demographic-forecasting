import pandas as pd
import numpy as np
import geopandas as gpd
import pickle
# import boto3

# TODO: import precalculated mus and hbar for SSP projections up to 2100

def load_dataframe_with_hbar_from_csv(filename):
    """
    Normally, importing a dataframe from a csv works well, but not if that dataframe
    contains a vector, which in the case of DFFT is a vector representing a discretization
    of the H function. So, this will load such a dataframe as long as the vector is surrounded
    by quotation marks
    :param filename:
    :return: df: panda dataframe
    """
    df = pd.read_csv(filename)
    # Convert headaches from string format into numpy array format
    def convert_string_arr_to_nparray(string_arr):
        # H = np.array(H_string[1:-1].split(', ')).astype(np.float)
        # print(H_string)
        string_arr = string_arr.replace('  ', ' ')
        string_arr = string_arr.replace('  ', ' ')
        string_arr = string_arr.replace('  ', ' ')
        string_arr = string_arr.replace('\n', '')
        string_arr = string_arr.replace(' ]', ']')
        string_arr = string_arr.replace('[ ', '[')
        if "," in string_arr:
            string_arr = string_arr[1:-1].split(', ')
        else:
            string_arr = string_arr[1:-1].split(' ')
        arr = np.array(string_arr).astype(np.float)
        return arr
    df["ns"] = [convert_string_arr_to_nparray(row.ns) for index, row in df.iterrows()]
    df["hbar"] = [convert_string_arr_to_nparray(row.hbar) for index, row in df.iterrows()]
    return df


def add_columns_to_blockgroup_df(df_blockgroup):
    # df_blockgroup = pd.read_pickle("blockgroup_longitudinal_data")
    df_blockgroup['fraction_black_90'] = df_blockgroup.black_90.div(df_blockgroup.total_90)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_black_90']), 'fraction_black_90'] = np.nan
    df_blockgroup['fraction_black_00'] = df_blockgroup.black_00.div(df_blockgroup.total_00)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_black_00']), 'fraction_black_00'] = np.nan
    df_blockgroup['fraction_black_10'] = df_blockgroup.black_10.div(df_blockgroup.total_10)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_black_10']), 'fraction_black_10'] = np.nan
    df_blockgroup['delta_fraction_black_00_10'] = df_blockgroup['fraction_black_10'] - df_blockgroup[
        'fraction_black_00']
    df_blockgroup['delta_fraction_black_90_10'] = df_blockgroup['fraction_black_10'] - df_blockgroup[
        'fraction_black_90']
    df_blockgroup['delta_fraction_black_90_00'] = df_blockgroup['fraction_black_00'] - df_blockgroup[
        'fraction_black_90']

    df_blockgroup['fraction_hispanic_90'] = df_blockgroup.hispanic_90.div(df_blockgroup.total_90)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_hispanic_90']), 'fraction_hispanic_90'] = np.nan
    df_blockgroup['fraction_hispanic_00'] = df_blockgroup.hispanic_00.div(df_blockgroup.total_00)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_hispanic_00']), 'fraction_hispanic_00'] = np.nan
    df_blockgroup['fraction_hispanic_10'] = df_blockgroup.hispanic_10.div(df_blockgroup.total_10)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_hispanic_10']), 'fraction_hispanic_10'] = np.nan
    df_blockgroup['delta_fraction_hispanic_00_10'] = df_blockgroup['fraction_hispanic_10'] - df_blockgroup[
        'fraction_hispanic_00']
    df_blockgroup['delta_fraction_hispanic_90_10'] = df_blockgroup['fraction_hispanic_10'] - df_blockgroup[
        'fraction_hispanic_90']
    df_blockgroup['delta_fraction_hispanic_90_00'] = df_blockgroup['fraction_hispanic_00'] - df_blockgroup[
        'fraction_hispanic_90']

    df_blockgroup['fraction_white_90'] = df_blockgroup.white_00.div(df_blockgroup.total_90)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_white_90']), 'fraction_white_90'] = np.nan
    df_blockgroup['fraction_white_00'] = df_blockgroup.white_00.div(df_blockgroup.total_00)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_white_00']), 'fraction_white_00'] = np.nan
    df_blockgroup['fraction_white_10'] = df_blockgroup.white_10.div(df_blockgroup.total_10)
    df_blockgroup.loc[~np.isfinite(df_blockgroup['fraction_white_10']), 'fraction_white_10'] = np.nan
    df_blockgroup['delta_fraction_white_00_10'] = df_blockgroup['fraction_white_10'] - df_blockgroup[
        'fraction_white_00']
    df_blockgroup['delta_fraction_white_90_10'] = df_blockgroup['fraction_white_10'] - df_blockgroup[
        'fraction_white_90']
    df_blockgroup['delta_fraction_white_90_00'] = df_blockgroup['fraction_white_00'] - df_blockgroup[
        'fraction_white_90']

    # df_blockgroup.to_pickle("./blockgroup_longitudinal_data_fractions.pickle")

    return df_blockgroup


def load_dataframe_with_H_from_csv(filename):
    """
    Normally, importing a dataframe from a csv works well, but not if that dataframe
    contains a vector, which in the case of DFFT is a vector representing a discretization
    of the H function. So, this will load such a dataframe as long as the vector is surrounded
    by quotation marks
    :param filename:
    :return: df: panda dataframe
    """
    df = pd.read_csv(filename)
    def convert_H_string_to_nparray(H_string):
        H = np.array(H_string[1:-1].split(', ')).astype(np.float)
        return H
    df["H"] = [convert_H_string_to_nparray(row.H) for index, row in df.iterrows()]
    return df


def load_county_GIS(filename = './data/cb_2018_us_county_20m.pickle'):
    with open(filename, 'rb') as f:
        geojson = pickle.load(f)
    geo_df_county = gpd.GeoDataFrame.from_features(geojson)
    return geo_df_county


def load_blockgroup_GIS(geoid_county,foldername="./data/block_group_shapefiles/"):
    if geoid_county<10000:
        filename = foldername + "0" + str(geoid_county) + '.pickle'
    else:
        filename = foldername + str(geoid_county) + '.pickle'

    # Read in geojson dictionary stored in a pickle file
    with open(filename, 'rb') as f:
        geojson_blockgroup = pickle.load(f)
    geo_df_blockgroup = gpd.GeoDataFrame.from_features(geojson_blockgroup)
    geo_df_blockgroup['GEOID10'] = pd.to_numeric(geo_df_blockgroup['GEOID10'])
    return geo_df_blockgroup, geojson_blockgroup


def load_blockgroup_GIS_S3(geoid_county,bucketname="neighborhood-forecasting-data"):
    # Loads the file from Amazon Bucket
    s3 = boto3.client('s3')
    s3.download_file('BUCKET_NAME', 'OBJECT_NAME', 'FILE_NAME')
    filename = geoid_county + '.pickle'
    # Read in geojson dictionary stored in a pickle file
    with open(filename, 'rb') as f:
        geojson_blockgroup = pickle.load(f)
    geo_df_blockgroup = gpd.GeoDataFrame.from_features(geojson_blockgroup)
    return geo_df_blockgroup, geojson_blockgroup



def convert_county_shp_to_geojson(county_filename='./data/cb_2018_us_county_20m.shp'):
    geo_county = gpd.read_file(county_filename)
    # Only extract the necessary information
    geojson_county = json.loads(geo_county.iloc[:, [0, 4, 5, 9]].to_json())  # Convert to String like object.
    # Now it is a dictionary. Now make the keys for each county to the GEOID
    for i in range(len(geojson_county['features'])):
        geojson_county['features'][i]['id'] = \
            geojson_county['features'][i]['properties']['GEOID']
    # Write to pickle (Preferred method)
    with open('cb_2018_us_county_20m.pickle', 'wb') as f:
        pickle.dump(geojson_county, f)

    return geojson_county


def convert_blockgroup_shps_to_geojson(county_filename='cb_2018_us_county_20m.shp'):
    geo_blockgroup = gpd.read_file(
        '/home/yunuskink/Desktop/nhgis0016_shape/US_blck_grp_2010.shp')
    geo_blockgroup = geo_blockgroup.to_crs(epsg=4326)
    # Only extract the necessary information
    shapefiles_for_2010 = 1
    if shapefiles_for_2010:
        geo_blockgroup["FIPSSTCO"] = geo_blockgroup["STATEFP10"] + geo_blockgroup["COUNTYFP10"]

    county_geoids = geo_blockgroup.FIPSSTCO.unique()
    for geoid in county_geoids:
        print(geoid)
        geojson_blockgroup = json.loads(geo_blockgroup.loc[geo_blockgroup["FIPSSTCO"] == geoid].to_json())
        for i in range(len(geojson_blockgroup['features'])):
            if shapefiles_for_2010:
                geojson_blockgroup['features'][i]['id'] = geojson_blockgroup['features'][i]['properties']['GEOID10']
            else:
                geojson_blockgroup['features'][i]['id'] = geojson_blockgroup['features'][i]['properties']['STFID']
        with open(
                '/home/yunuskink/Desktop/nhgis0016_shape/county_subsets' + geoid + '.pickle',
                'wb') as f:
            pickle.dump(geojson_blockgroup, f)
    # geojson_blockgroup = json.loads(geo_blockgroup.iloc[:, [0, 4, 5, 9]].to_json())  # Convert to String like object.
    # Now it is a dictionary. Now make the keys for each county to the GEOID
    return

def load_global_f_V(filename = "./parameters/headaches_yunus.csv"):
    global_fit_df = load_dataframe_with_H_from_csv(filename)
    races = ["black","hispanic","white"]
    decades = [1990, 2000, 2010]
    f_dict = {}
    ns = np.linspace(0,1,len(global_fit_df.iloc[0]["H"]))
    rows = []
    for race in races:
        for decade in decades:
            decade_str = str(decade)
            f = global_fit_df[(global_fit_df["decade"] == decade) &
                              (global_fit_df["race"] == race) &
                              (global_fit_df["county"] == 1001)]["H"].iloc[0]
            f = f - f[0]
            f = f - ns*f[-1]
            f_dict[race+"_"+decade_str[2:]] = f
            df_subset = global_fit_df[(global_fit_df["decade"] == decade) &
                                      (global_fit_df["race"] == race)]
            for row in df_subset.itertuples(index=False):
                H = row.H
                V = H[-1] - f[0] - f[-1]
                V_row = {"decade": decade, "race": race, "county": row.county, "V": V}
                rows.append(V_row)
    V_df = pd.DataFrame(rows)
    return f_dict, V_df


