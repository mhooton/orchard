import pandas as pd
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import griddata
from scipy import signal, ndimage
import matplotlib.pyplot as plt
import io
import requests
from astropy.time import Time
import json
# from joblib import Memory

# from julian.julian import to_jd, from_jd

import time, sys

# memory = Memory(location="./cache/", verbose=0)
# memory1 = Memory(location="./cache/", verbose=0)
# memory2 = Memory(location="./cache/", verbose=0)

plt.style.use('seaborn-v0_8-notebook')

width, height = plt.figaspect(5. / 8.)


def interpolate_dfs(wavelengths, *data):
    '''
    Interpolates panda dataframes onto an index, of same index type (e.g. wavelength in microns)

    Parameters
    ----------
    wavelength: 1d array which data is to be interpolated onto
    data:       Pandas dataframes

    Returns
    -------
    df: Interpolated dataframe

    '''
    df = pd.DataFrame({'tmp': wavelengths}, index=wavelengths)
    for dat in data:
        df = pd.concat([df, dat], axis=1)
    df = df.interpolate('index').reindex(wavelengths)
    df = df.drop('tmp', 1)
    # df = df.fillna(0)
    return df


def generateBase(sResponse, basedir=None):
    '''
    Generates the water vapour grid base for Paranal, Chile. Takes a few minutes.

    Generates a base grid for:
    airmass: 1 - 3
    pwv: 0.05 - 30 mm
    Teff: 2000 - 36500 K

    See arrays for base resolutions

    Ref: ...

    Parameters
    ----------
    sResponse:  csv file with two (unlabelled) columns, wavelength (in microns), system spectral response curves of telescope + filter + camera (as fraction).

    Returns
    -------
    coords, data: coordinates and data of base grid generated.

    '''
    if basedir is None:
        # Try to find it relative to current location
        basedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    grid_ingredients_path = os.path.join(basedir, "pwv", "datafiles", "gridIngredients.pkl")

    if not os.path.exists(grid_ingredients_path):
        raise FileNotFoundError(f"gridIngredients.pkl not found at {grid_ingredients_path}")

    gridIngredients = pd.read_pickle(grid_ingredients_path)
    rsr = pd.read_csv(sResponse, header=None, index_col=0)

    pwv_values = np.array([0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 7.5, 10.0, 20.0, 30.0])
    airmass_values = np.array(
        [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0])
    temperature_values = np.array(
        [2000, 2100, 2250, 2320, 2400, 2440, 2500, 2600, 2650, 2710, 2850, 3000, 3030, 3100, 3200, 3250, 3410, 3500,
         3550, 3650, 3700, 3800, 3870, 3940, 4000, 4070, 4190, 4230, 4330, 4410, 4540, 4600, 4700, 4830, 4990, 5040,
         5140, 5170, 5240, 5280, 5340, 5490, 5530, 5590, 5660, 5680, 5720, 5770, 5880, 5920, 6000, 6060, 6170, 6240,
         6340, 6510, 6640, 6720, 6810, 7030, 7220, 7440, 7500, 7800, 8000, 8080, 8270, 8550, 8840, 9200, 9700, 10400,
         10700, 12500, 14000, 14500, 15700, 16700, 17000, 18500, 20600, 24500, 26000, 29000, 31500, 32000, 32500, 33000,
         34500, 35000, 36500])

    wavelengths = np.arange(0.5, 2, 0.0001)

    gridSauce = interpolate_dfs(wavelengths, rsr, gridIngredients)
    gridSauce = gridSauce[(gridSauce[1] > 0)]
    atm_grid = []
    for i, pwv in enumerate(pwv_values):
        # print(pwv)
        update_progress(i / (len(pwv_values) - 1))
        for airmass in airmass_values:
            for temperature in temperature_values:
                atmosphere_trans = gridSauce[str(pwv) + '_' + str(airmass)]
                simStar = gridSauce[str(temperature) + 'K']
                response = simps(gridSauce[1] * atmosphere_trans * simStar / max(simStar), gridSauce.index)

                atm_grid.append((pwv, airmass, temperature, response))

    data = np.array([x[3] for x in atm_grid])
    data = data.reshape((len(pwv_values), len(airmass_values), len(temperature_values)))

    coords = np.zeros((len(pwv_values), len(airmass_values), len(temperature_values), 3))
    coords[..., 0] = pwv_values.reshape((len(pwv_values), 1, 1))
    coords[..., 1] = airmass_values.reshape((1, len(airmass_values), 1))
    coords[..., 2] = temperature_values.reshape((1, 1, len(temperature_values)))

    # np.save('./water_grid_coords_20191023.npy', coords)
    # np.save('./water_grid_I+z_data_20191023.npy', data)

    return coords, data


# @memory.cache
def interp(coords, data, pwv, airmass, Teff):
    '''
    Interpolates between water grid base points (waterGrid.generateBase(...)), using a cubic method.

    Parameters
    ----------
    coords, data:   coordinates and data of base grid generated.
    pwv:            precipitable water vapour value at zenith
    airmass:        airmass of target/comparison star
    Teff:           effective temperature of target/comparison star

    Returns
    -------
    interp: interpolated value of grid.

    '''
    method = 'cubic'
    Teffs = coords[..., 2][0, 0]
    Teff_lower = np.max(Teffs[Teffs <= Teff])
    Teff_upper = np.min(Teffs[Teffs >= Teff])

    if Teff_lower == Teff_upper:
        x = coords[..., 0][coords[..., 2] == Teff]  # pwv
        y = coords[..., 1][coords[..., 2] == Teff]  # airmass
        z = data[coords[..., 2] == Teff]  # effect

        interp = griddata((x, y), z, (pwv, airmass), method=method)  # interpolated value
    else:
        x_lower = coords[..., 0][coords[..., 2] == Teff_lower]  # pwv
        y_lower = coords[..., 1][coords[..., 2] == Teff_lower]  # airmass
        z_lower = data[coords[..., 2] == Teff_lower]  # effect
        interp_lower = griddata((x_lower, y_lower), z_lower, (pwv, airmass),
                                method=method)  # interpolated value lower Teff

        x_upper = coords[..., 0][coords[..., 2] == Teff_upper]  # pwv
        y_upper = coords[..., 1][coords[..., 2] == Teff_upper]  # airmass
        z_upper = data[coords[..., 2] == Teff_upper]  # effect
        interp_upper = griddata((x_upper, y_upper), z_upper, (pwv, airmass),
                                method=method)  # interpolated value upper Teff

        w_lower = (Teff_upper - Teff) / (Teff_upper - Teff_lower)  # lower weight
        w_upper = (Teff - Teff_lower) / (Teff_upper - Teff_lower)  # upper weight

        interp = w_lower * interp_lower + w_upper * interp_upper  # final interpolated value

    return interp


def deltaFluxPlot(coords, data, name=''):
    '''
    Comment and add more flexibility?

    '''
    cm = plt.cm.get_cmap('jet_r')
    airmass = 1
    pwv_values = np.arange(0.05, 10.55, 0.5)
    temperatures = ticks = [6000, 5500, 5000, 4500, 4000, 3500, 3000, 2500, 2000]
    iticks = ['G0', 'G7', 'K2', 'K5', 'K8', 'M3', 'M5', 'M8', 'L2']

    eff_df = pd.DataFrame({'pwv': pwv_values})

    fig, ax = plt.subplots(figsize=(width, height), dpi=150)

    for k, temperature in enumerate(temperatures):
        val = interp(coords, data, min(pwv_values), airmass, temperature)
        eff = []
        for pwv in pwv_values:
            val2 = interp(coords, data, pwv, airmass, temperature)
            eff.append(100 * (val2 - val) / val)

        eff_df[temperature] = eff

    for k, temperature in enumerate(temperatures):
        sc = plt.scatter(eff_df['pwv'], eff_df[temperature], c=np.ones(len(eff_df['pwv'])) * temperature, cmap=cm,
                         alpha=0.8, vmin=min(temperatures), vmax=max(temperatures), marker='None', label='_nolegend_')

        ax.plot(eff_df['pwv'], eff_df[temperature], color=sc.to_rgba(temperature),
                label=iticks[k] + ", " + str(temperature) + ' K')

    cbar = plt.colorbar(sc);
    cbar2 = plt.colorbar(sc)

    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(iticks)

    cbar2.set_ticks(ticks)
    cbar2.ax.set_yticklabels(ticks)

    cbar.set_label(r'Spectral Type')
    cbar2.set_label(r'T$_{\rm eff}$ [K]')

    ax.legend(framealpha=0.2)

    ax.set_xticks(np.arange(0, 11, 1))
    ax.set_ylim(-10.5, 0.5)
    ax.set_xlabel('PWV [mm]')
    ax.set_ylabel(r'$\Delta$flux [\%]')
    #ax.set_title('Name: \'%s\'' % (name))
    ax.minorticks_on()

    return fig, ax

# @memory2.cache
def getLHATPROdata(dateS, dateE, peakR=False, format='jd'):
    data = {
        'wdbo': 'tsv/display',
        'max_rows_returned': '50000',
        'start_date': dateS + '..' + dateE,
        'tab_integration': 'on',
        'integration': '',
        'tab_irt0': 'on',
        'irt0': '',
        'tab_lwp0': 'on',
        'lwp0': '',
        'tab_pwv0': 'on',
        'pwv0': '',
        'order': 'start_date'
    }
    print(("Getting LHATPRO data for:", dateS, dateE))
    response = requests.post('http://archive.eso.org/wdb/wdb/asm/lhatpro_paranal/query', data=data)

    df = pd.read_csv(io.StringIO(response.text), delimiter='\t')

    if peakR:
        indexes = df[df['IR temperature [Celsius]'] == '                              '].index
        # Apply filtering to a COPY of the dataframe
        df_filtered = df.copy()
        df_filtered = df_filtered.drop(indexes[:-2] + 1)
        indexes = df_filtered[df_filtered['IR temperature [Celsius]'] == '                              '].index
        df_filtered = df_filtered.drop(indexes)

        # Extract timestamps and PWV from the FILTERED dataframe
        try:
            time = pd.to_datetime(df_filtered['Date time'].str.strip(), format='%Y-%m-%dT%H:%M:%S')
        except Exception as e:
            return None

        # Extract PWV values from the FILTERED dataframe
        pwv_values = df_filtered['Precipitable Water Vapour [mm]'].str.strip()
        pwv_values = pd.to_numeric(pwv_values, errors='coerce')
        pwv_array = pwv_values.values

        pwvData = pd.DataFrame({'pwv': pwv_array}, index=time)

        if (format == 'jd'):
            pwvData.index = Time(pwvData.index.values).jd
            pwvData_orig.index = Time(pwvData_orig.index.values).jd

        pwvData = peak_removal(pwvData)

        # Return both original and processed data
        pwvData.pwv_original = pwvData_orig.pwv.values
        pwvData.time_original = pwvData_orig.index.values

    else:
        # Parse timestamps first and check for failures
        try:
            time = pd.to_datetime(df['Date time'].str.strip(), format='%Y-%m-%dT%H:%M:%S', errors='coerce')
        except Exception as e:
            return None

        # Clean and convert PWV values
        pwv_values = df['Precipitable Water Vapour [mm]'].str.strip()
        pwv_values = pd.to_numeric(pwv_values, errors='coerce')

        # Create DataFrame and remove rows where EITHER timestamp OR PWV is invalid
        valid_mask = ~(time.isna() | pwv_values.isna())
        if valid_mask.sum() == 0:
            print("ERROR: No rows with both valid timestamp and PWV")

            return None

        # Filter to valid rows only
        time_clean = time[valid_mask]
        pwv_clean = pwv_values[valid_mask]

        # Create DataFrame with clean data
        pwvData = pd.DataFrame({'pwv': pwv_clean.values}, index=time_clean)

        if (format == 'jd'):
            try:
                pwvData.index = Time(pwvData.index.values).jd
            except Exception as e:
                return None

    return pwvData


# @memory1.cache
def effect(coords, data, jd, airmass, temp, sigma=0):
    '''
    Worried about integration time, when does a PWV value match an observation time?

    '''
    dateS = Time(np.floor(min(jd)) - 0.5, format='jd').datetime.strftime("%Y-%m-%d")
    dateE = Time(np.ceil(max(jd)) + 0.5, format='jd').datetime.strftime("%Y-%m-%d")

    pwvData = getLHATPROdata(dateS, dateE)

    df = interpolate_dfs(jd, pwvData)
    if sigma != 0:
        df.pwv = ndimage.gaussian_filter1d(df.pwv, sigma)

    target = interp(coords, data, df.pwv, airmass, temp)

    targetdf = pd.DataFrame({'normFlux': target / (np.mean(target)), 'pwv': df.pwv, 'airmass': airmass.values},
                            index=jd)
    return targetdf


def getATM(pwv, airmass):
    gridIngredients = pd.read_pickle('./datafiles/gridIngredients.pkl')
    return gridIngredients[str(pwv) + '_' + str(airmass)]


def generateSR(efficiencyFile, filterFile, SRFile):
    wavelengths = np.arange(0.5, 2, 0.0001)

    eff = pd.read_csv(efficiencyFile, header=None)
    filt = pd.read_csv(filterFile, header=None)

    effDF = pd.DataFrame({'eff': eff[1].values}, index=eff[0])
    filtDF = pd.DataFrame({'filt': filt[1].values}, index=filt[0])
    df = interpolate_dfs(wavelengths, effDF, filtDF)

    dfSR = df['eff'] * df['filt']

    dfSR = dfSR[dfSR > 0]

    dfSR.to_csv(SRFile, header=False)

    dfSR.plot()

    print((SRFile + " has been saved."))


def update_progress(progress):
    """Simple progress indicator without IPython dependency"""
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))

    # Simple progress bar that works in any terminal
    text = "\rProgress: [{0}] {1:.1f}%".format("#" * block + "-" * (bar_length - block), progress * 100)
    print(text, end='', flush=True)

    if progress >= 1:
        print()  # New line when complete

def peak_removal(pwvData):
    peak_widths = np.arange(1, 3)
    peak_indices = signal.find_peaks_cwt(pwvData.pwv, peak_widths)

    ind_array = []
    for i in peak_indices:
        ind_array.append(i + 1)
        ind_array.append(i - 1)
        ind_array.append(i)

    # NEW:
    jd_list = []
    pwvR = []
    for i, val in enumerate(pwvData.index):
        if i not in ind_array:
            jd_list.append(pwvData.index[i])
            pwvR.append(pwvData.pwv.iloc[i])

    return pd.DataFrame({'pwv': pwvR}, index=jd_list)
