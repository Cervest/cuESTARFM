import numpy as np
import glob
from datetime import datetime as dt
import xml.etree.ElementTree as ET
import matplotlib.dates as dates
import datetime
import subprocess as sp
import os


def doy2date(year, doy):
    """
    Convert day of year to date if year is known
    """
    date_0 = datetime
    jul_0 = dates.num2julian( dates.date2num(date_0.date(year,1,1)) )
    date = dates.num2date(dates.julian2num(jul_0+doy-1)) #jul2date(jul_0+doy-1)

    return date


def get_sort_l8(lst_xml, band_position):

    date_arr = []
    band_arr = []
    for fname in lst_xml:
        tree = ET.parse(fname)
        root = tree.getroot()
        date_str = root[0][3].text
        date_arr = np.append(date_arr, dt.strptime(date_str, '%Y-%m-%d').date())
        path = os.path.dirname(fname)
        fname_band = path + '/' + root[1][band_position][2].text
        print(fname_band)
        band_arr = np.append(band_arr, fname_band)
    date_ind = np.argsort(date_arr)
    date_arr = date_arr[date_ind]
    lst_sort = lst_xml[date_ind]
    band_arr = band_arr[date_ind]

    return lst_sort, date_arr, band_arr


def do_gdalwarp(file_shp, fname_in, fname_out):

    cmd = f'gdalwarp -t_srs EPSG:4326  -crop_to_cutline -cutline  {file_shp} {fname_in} {fname_out}'
    sp.check_output(cmd, shell=True)


def get_lst_l8(file_l8, dir_out, file_shp, band, date_start, date_end):

    lst = np.array(glob.glob(file_l8))

    lst_sort, date_l8_all, band_arr = get_sort_l8(lst, band)
    # print(band_arr),

    l_arr = []
    date_l8 = []
    for i, fname in enumerate(band_arr):
        if np.logical_and(date_l8_all[i] > date_start, date_l8_all[i] < date_end):
            date_l8 = np.append(date_l8, date_l8_all[i])
            fname_out = date_l8[-1].strftime(f'{dir_out}l_%Y_%m_%d_b{band - 1}.TIF')
            l_arr = np.append(l_arr, fname_out)
            do_gdalwarp(file_shp, fname, fname_out)

    return l_arr, date_l8


def get_lst_modis(file_modis, dir_out, file_shp, band, date_start, date_end):

    lst_mod = glob.glob(file_modis % band)
    m_arr = []
    date_mod = []
    for i, fname in enumerate(lst_mod):
        date = dt.strptime(fname.split('/')[-2], '%Y%j').date()
        if np.logical_and(date > date_start, date < date_end):
            print(fname)
            date_mod = np.append(date_mod, date)
            fname_out = date.strftime(f'{dir_out}m_%Y_%m_%d_b{band}.TIF')
            do_gdalwarp(file_shp, fname, fname_out)
            m_arr = np.append(m_arr, fname_out)

    return m_arr, date_mod


if __name__ == '__main__':

    bands_landsat = np.array([3, 4, 5, 6, 7, 8])
    bands_modis =   np.array([3, 4, 1, 2, 6, 7])

    file_shp = '../data/vector/fr_rect_small/fr_rect_small.shp'
    file_l8 = '/data/landsat-8/france/LC08*2018*T1-*/*xml'
    file_modis = '/data/raw_data/modis/mcd43a4_aws/18/04/2018*/*.A2018*.*_B%02d.TIF'
    dir_out = '../data/intermediate/'

    date_start = datetime.date(2018, 4, 20)
    date_end = datetime.date(2018, 9, 13)

    sp.check_output(f'rm {dir_out}*', shell=True)

    # l_arr, date_l8 = get_lst_l8(file_l8, dir_out, bands_landsat, date_start, date_end)
    # m_arr, date_mod = get_lst_modis(file_modis, dir_out, bands_modis, date_start, date_end)

    # lst_mod = []
    # for fname in lst_mod0:
    #     date = dt.strptime(fname.split('/')[-2], '%Y%j').date()
    #     if np.logical_and(date>datetime.date(2018,7,17), date<datetime.date(2018,8,15)):
    #         lst_mod.append(fname)

    # lst = np.array(glob.glob(file_l8))

    for band_i in range(bands_landsat.shape[0]):

        l_arr, date_l8 = get_lst_l8(file_l8, dir_out, file_shp, bands_landsat[band_i], date_start, date_end)
        m_arr, date_mod = get_lst_modis(file_modis, dir_out, file_shp, bands_modis[band_i], date_start, date_end)

        # Write strings in-pair MODIS and Landsat. I.e. write pairs of hi-low resolution which
        # will be used for prediction
        in_pair_modis = ''
        in_pair_landsat = ''
        ind_mod_pair = []
        for i, d in enumerate(date_l8):
            d_day = [dd.days for dd in (date_mod - d)]
            ind_mod_pair = np.append(ind_mod_pair, np.argmin(np.abs(d_day))).astype(int)
            # print(l_arr[i], m_arr[ind])
            in_pair_landsat = in_pair_landsat + l_arr[i] + ' '
            in_pair_modis = in_pair_modis + m_arr[ind_mod_pair[-1]] + ' '
        ind_mod = np.arange(m_arr.shape[0])
        ind_mod = np.delete(ind_mod, ind_mod_pair).astype(int)

        in_pday_modis = ''
        for i in range(0, ind_mod.shape[0]):
            in_pday_modis = in_pday_modis + m_arr[ind_mod[i]] + ' '

        param_begin = 'ESTARFM_PARAMETER_START\n\n'
        in_pair_modis = in_pair_modis.replace('intermediate', 'results').replace('m_', 'M_') + '\n\n'
        in_pair_landsat = in_pair_landsat.replace('intermediate', 'results').replace('l_', 'L_') + '\n\n'
        in_pday_modis = in_pday_modis.replace('intermediate', 'results').replace('m_', 'M_') + '\n\n'
        in_pday_landsat = in_pday_modis.replace('M_', 'LM_').replace(f'b{bands_modis[band_i]}', f'b{bands_landsat[band_i]-1}') + '\n\n'
        num_in_pairs =  f'NUM_IN_PAIRS = {l_arr.shape[0]}\n\n'

        param_end = 'The_width_of_searching_window = 51\n' +\
                    'Assumed_number_of_classifications = 6\n' +\
                    'sensor_uncertain = 0.0028\n' +\
                    'NODATA = 0\n' +\
                    'G_Type = GTIff\n\n' +\
                    'ESTARFM_PARAMETER_END\n'

        f = open(f'../Codes/param_max_b{bands_landsat[band_i]-1}.txt', 'w')
        f.write(param_begin +
                num_in_pairs +
                f'IN_PAIR_MODIS_FNAME = {in_pair_modis}' +
                f'IN_PAIR_LANDSAT_FNAME =  {in_pair_landsat}' +
                f'IN_PDAY_MODIS_FNAME = {in_pday_modis}' +
                f'OUT_PDAY_LANDSAT_FNAME = {in_pday_landsat}' +
                param_end)
        f.close()

        l_str = ' '.join(l_arr)
        m_str = ' '.join(m_arr)
        f_merged = f'{dir_out}merged.tif'
        cmd = f'gdal_merge.py -of gtiff -separate -o {f_merged} {l_str} {m_str}'

        print(cmd)
        sp.check_output(cmd, shell=True)

        n = len(l_arr) + len(m_arr)
        all_arr = np.concatenate((l_arr, m_arr))
        for i in range(n):
            # fname_out = date.strftime('../data/results/m_%Y_%m_%d.TIF')
            fname_out = all_arr[i].replace('intermediate', 'results').replace('l_', 'L_').replace('m_', 'M_')
            sp.check_output(f'gdal_translate -b {i+1} {f_merged} {fname_out}', shell=True)
