# FM 1 min raw data processing code
# Created by: Diane Wang (tdwang@ucdavis.edu)
# Version 0
# Version Date: 6/30/2025
# Updated Date: 7/2/2025

# About:
# This code is designed to process the raw data from the FM field campaign. Both the LICOR and the CR1000x systems have 1-min slow measurements.
# The code first unzips the LICOR raw data and then reads data from both systems and stores into one file for future analysis.

# Steps for use:
#
# 1. Modify paths for both raw data and output data
# 2. Unzip LICOR system data
# 3. Use pandas to read data and store into a new dataframe; modify the stored variables as needed.
# Note: Only NaN values are filtered. Future filters can also be applied.


import os
import pandas as pd
from zipfile import ZipFile
import glob

# *****1. Modify paths below
c_path = os.getcwd()
if 'diane_wt' in c_path:
    box_path = '/Users/diane_wt/Library/CloudStorage/Box-Box'
elif 'admin-dianew68' in c_path:
    box_path = 'C:/Users/admin-dianew68/Box'

# cr1000xFM_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250527_tall/ascii/' # FM Main Tower data dir, need to change depending on processing dates
# cr1000xFM_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250828_tall/ascii/' # FM Main Tower data dir, need to change depending on processing dates
cr1000xFM_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20251110_tall/ascii/' # FM Main Tower data dir, need to change depending on processing dates
FM_processed_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/processed/FM_DOL_1min/' # output data dir

# licor_zipfile_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250527_DOL50m_Data/DOL DAqM/' # zip file folders for raw 1-min LICOR data
licor_zipfile_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250828_DOL50m_Data/2025-11-10/Daqm/'
# unzipfile_dir = box_path + '/Diane/French Meadows/data/DOL_unziped/20250401-20250418/' # path to save unzipped files for LICOR data, this folder can be your local or other box folder
# unzipfile_dir = box_path + '/Diane/French Meadows/data/DOL_unziped/20250101-20250331/' # path to save unzipped files for LICOR data, this folder can be your local or other box folder
unzipfile_dir = box_path + '/Diane/French Meadows/data/DOL_unziped/20250901-20251110/' # path to save unzipped files for LICOR data, this folder can be your local or other box folder

# *****2. Unzip LICOR date
# moni = [4]
# daylength = 18 # unzip data day length, example here is from 20250401 to 20250418, adds up to 18 days
# moni = [1, 2, 3]
# daylength = [31, 28, 31]
# moni = [11]
# daylength = [9]
# for mi in range(len(moni)):
#     for di in range(daylength[mi]): # extract one more because the 00:00:00 record is written in the next daily file
#         zipfile = licor_zipfile_dir + f'2025-{moni[mi]:02}-{1+di:02}-000000-daqm.zip'  # ****change here if you want to start from another day
#         z = ZipFile(zipfile)
#         z.extractall(unzipfile_dir)

log_file_path = glob.glob(os.path.join(unzipfile_dir, f'2025-*-000000-daqm.log'))
log_file_path.sort()
filelength = len(log_file_path)
# *****3. Combine data tables and filter NaN values
for fi in range(filelength-1): # 1
    # the following code read two daily files from LICOR, because LICOR data starts from 00:00:00 of the day but CR1000x data starts from 00:00:01
    licor_df = pd.read_csv(log_file_path[fi], sep='\s+', engine='python', header=0,
                           skiprows=[1, 2])  # skip first row of unit, second row of the 00:00:00 record
    licor_df_nextfirstline = pd.read_csv(log_file_path[fi+1], sep='\s+', engine='python', header=0,
                                         skiprows=[1], na_values=['-9999'])  # skip the first row of unit, keep the 00:00:00 record for the next day
    licor_df = pd.concat([licor_df, licor_df_nextfirstline.iloc[[0]]], ignore_index=True)  # create a new table to be consistent with CR1000x data table
    # create indices using timestamps and drop duplicates
    licor_df['TIMESTAMP'] = pd.to_datetime(licor_df['DATE'] + ' ' + licor_df['TIME'])
    licor_df.set_index('TIMESTAMP', inplace=True)
    licor_df = licor_df[~licor_df.index.duplicated(keep='first')]

    # cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{516+fi}.dat" # 4/1/2025-4/18/2025
    # cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{426 + fi}.dat"  # 1/1/2025-3/31/2025
    # cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{638 + fi}.dat"  # 7/24/2025-8/28/2025
    if 681 + fi<= 718:
        cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{681 + fi}.dat"  # 9/1/2025-10/08/2025
    else:
        cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{681 + fi + 1}.dat"  # 10/09/2025-11/09/2025
    cr1000x_df = pd.read_csv(cr1000x_file_path, skiprows=[0, 2, 3], index_col=[0], na_values=['NaN', 'NAN'],
                             parse_dates=True)  # read data

    # cr1000x_file_path = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{636 + fi}.dat"      # ****** this code only used for 7/23/2025
    # cr1000x_file_path2 = f"{cr1000xFM_dir}TOA5_6653slow_avg_data{636 + fi + 1}.dat" # ****** this code only used for 7/23/2025
    # cr1000x_files = [cr1000x_file_path, cr1000x_file_path2]
    # cr1000x_df = pd.concat(
    #     [pd.read_csv(f, skiprows=[0, 2, 3], index_col=[0], na_values=['NaN', 'NAN'], parse_dates=True) for f in
    #      cr1000x_files
    #      ], axis=0)                                                                 # ****** this code only used for 7/23/2025

    cr1000x_df = cr1000x_df[~cr1000x_df.index.duplicated(keep='first')]
    cr1000x_df = cr1000x_df.reindex(licor_df.index)

    df = pd.DataFrame({})
    df.insert(loc=0, column='year', value=[int(stamp.year) for stamp in licor_df.index])
    df.insert(loc=1, column='day', value=[int(stamp.day_of_year) for stamp in licor_df.index])
    df.insert(loc=2, column='HM', value=[int(f"{stamp.hour:02d}{stamp.minute:02d}") for stamp in licor_df.index])
    df.insert(loc=3, column='second', value=0)

    df.insert(loc=4, column='Temp_6.35', value=cr1000x_df['AirTC_ee181_1_Avg'].values)
    df.insert(loc=5, column='RH_6.35', value=cr1000x_df['RH_ee181_1_Avg'].values)
    df.insert(loc=6, column='Temp_15', value=licor_df['TA_1_1_1'].values)
    df.insert(loc=7, column='RH_15', value=licor_df['RH_1_1_1'].values)
    df.insert(loc=8, column='Temp_30', value=licor_df['TA_1_2_1'].values)
    df.insert(loc=9, column='RH_30', value=licor_df['RH_1_2_1'].values)
    df.insert(loc=10, column='Temp_51.5', value=licor_df['TA_1_3_1'].values)
    df.insert(loc=11, column='RH_51.5', value=licor_df['RH_1_3_1'].values)
    df.insert(loc=12, column='SWIN_4.22', value=cr1000x_df['SWTop_Avg'].values)
    df.insert(loc=13, column='SWOUT_4.22', value=cr1000x_df['SWBottom_Avg'].values)
    df.insert(loc=14, column='LWIN_4.22', value=cr1000x_df['LWTopC_Avg'].values)
    df.insert(loc=15, column='LWOUT_4.22', value=cr1000x_df['LWBottomC_Avg'].values)
    df.insert(loc=16, column='Rn_4.22', value=cr1000x_df['Rn_Avg'].values)
    df.insert(loc=17, column='ALBEDO_4.22', value=cr1000x_df['albedo_Avg'].values)
    df.insert(loc=18, column='SWIN_51.5', value=licor_df['SWIN_1_1_1'].values)
    df.insert(loc=19, column='SWOUT_51.5', value=licor_df['SWOUT_1_1_1'].values)
    df.insert(loc=20, column='LWIN_51.5', value=licor_df['LWIN_1_1_1'].values)
    df.insert(loc=21, column='LWOUT_51.5', value=licor_df['LWOUT_1_1_1'].values)
    df.insert(loc=22, column='Rn_51.5', value=licor_df['RN_1_1_1'].values)
    df.insert(loc=23, column='ALBEDO_51.5', value=licor_df['ALB_1_1_1'].values)

    df.apply(pd.to_numeric)
    datestr = f'{cr1000x_df.index[0].year}{cr1000x_df.index[0].month:02d}{cr1000x_df.index[0].day:02d}'
    df.to_csv(
        FM_processed_dir + f'FM_DOL_1min_{datestr}000000.txt',
        header=True, index=None, sep=',')
    del licor_df, cr1000x_df, df