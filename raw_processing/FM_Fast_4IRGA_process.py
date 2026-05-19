import os
import pandas as pd

c_path = os.getcwd()
if 'diane_wt' in c_path:
    box_path = '/Users/diane_wt/Library/CloudStorage/Box-Box'
elif 'admin-dianew68' in c_path:
    box_path = 'C:/Users/admin-dianew68/Box'

FM_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250718_tall/ascii/' # FM Main Tower data dir
FM_processed_dir = box_path + "/Lab Library/French Meadows/Summer2025/data/processed/FM_DOL_10Hz/" # processed data dir

# 597, 599; 575, 597; 600, 627
startid = 604 # **** change this index corresponds to a start date
endid = 607   # **** change this index corresponds to an end date
HZ = [10]

# last_row_previous = None  # initialize outside loop
for fi in range(startid, endid, 2):
    cr1000x_file1 = f"{FM_dir}TOA5_6653fast_20hz_data{fi}.dat"
    cr1000x_file2 = f"{FM_dir}TOA5_6653fast_20hz_data{fi + 1}.dat"
    # cr1000x_file3 = f"{FM_dir}TOA5_6653fast_20hz_data{fi + 2}.dat" # apply this code if need more files to combine a 48-hour data file
    # cr1000x_files = [cr1000x_file1, cr1000x_file2, cr1000x_file3]  # apply this code if need more files to combine a 48-hour data file
    cr1000x_files = [cr1000x_file1, cr1000x_file2]
    cr1000x_df = pd.concat([pd.read_csv(f, skiprows=[0, 2, 3], index_col=[0], na_values=['NaN', 'NAN'], parse_dates=True) for f in cr1000x_files
], axis=0) # read data
    cr1000x_df = cr1000x_df[~cr1000x_df.index.duplicated(keep='first')]
    timestrings = pd.to_datetime(cr1000x_df.index, format='mixed')
    # ensure your dataframe index is datetime64[ns] (just in case)
    cr1000x_df.index = pd.to_datetime(timestrings)

    # Create full time index covering 48 hours at 0.1 s step
    start_time = pd.to_datetime(timestrings[0].floor('D')) + pd.Timedelta(seconds=1/HZ[0])  # midnight of first day at the first dt time record
    end_time = start_time + pd.Timedelta(hours=48)
    full_time_index = pd.date_range(start=start_time, end=end_time, freq=f'{1000/HZ[0]}ms', inclusive='left')
    # Ensure both are pandas datetime64[ns] type
    full_time_index = pd.to_datetime(full_time_index, format="%Y-%m-%d %H:%M:%S.%f")
    # Reindex your dataframe to this full time index
    cr1000x_df_full = cr1000x_df.reindex(full_time_index)
    # Continue your processing pipeline using cr1000x_df_full instead of cr1000x_df

    # # Only from second loop onwards: fill first row with last_row_previous if missing
    # if last_row_previous is not None:
    #     if cr1000x_df_full.iloc[0].isna().all():
    #         cr1000x_df_full.iloc[0] = last_row_previous.values

    timestrings = cr1000x_df_full.index  # update timestrings to full time index

    print(f"Expected rows: {len(full_time_index)}, Actual rows before reindex: {len(cr1000x_df)}, and after reindex: {len(cr1000x_df_full)}")

    df = pd.DataFrame({})
    df.insert(loc=0, column='year', value=[int(stamp.year) for stamp in timestrings])
    df.insert(loc=1, column='day', value=[int(stamp.day_of_year) for stamp in timestrings])
    df.insert(loc=2, column='HM', value=[int(f"{stamp.hour:02d}{stamp.minute:02d}") for stamp in timestrings])
    df.insert(loc=3, column='second', value=[f"{stamp.second + stamp.microsecond/1000000:.1f}" for stamp in timestrings]) # *****Change here if data is in 20 Hz*****

    df.insert(loc=4, column='Ux_32.18', value=cr1000x_df_full['Ux_2'].values)
    df.insert(loc=5, column='Uy_32.18', value=cr1000x_df_full['Uy_2'].values)
    df.insert(loc=6, column='Uz_32.18', value=cr1000x_df_full['Uz_2'].values)
    df.insert(loc=7, column='T_Sonic_32.18', value=cr1000x_df_full['Ts_2'].values)
    df.insert(loc=8, column='diagnostic_32.18', value=cr1000x_df_full['diag_sonic_2'].values)
    df.insert(loc=9, column='Pressure_32.18', value=cr1000x_df_full['cell_press_2'].values)
    df.insert(loc=10, column='H2O_32.18', value=cr1000x_df_full['H2O_2'].values)
    df.insert(loc=11, column='H2Osig_32.18', value=cr1000x_df_full['H2O_sig_strgth_2'].values)
    df.insert(loc=12, column='CO2_32.18', value=cr1000x_df_full['CO2_2'].values)
    df.insert(loc=13, column='CO2sig_32.18', value=cr1000x_df_full['CO2_sig_strgth_2'].values)
    df.insert(loc=14, column='gas_diag_32.18', value=cr1000x_df_full['diag_irga_2'].values)

    df.insert(loc=15, column='Ux_13.94', value=cr1000x_df_full['Ux_3'].values)
    df.insert(loc=16, column='Uy_13.94', value=cr1000x_df_full['Uy_3'].values)
    df.insert(loc=17, column='Uz_13.94', value=cr1000x_df_full['Uz_3'].values)
    df.insert(loc=18, column='T_Sonic_13.94', value=cr1000x_df_full['Ts_3'].values)
    df.insert(loc=19, column='diagnostic_13.94', value=cr1000x_df_full['diag_sonic_3'].values)
    df.insert(loc=20, column='Pressure_13.94', value=cr1000x_df_full['cell_press_3'].values)
    df.insert(loc=21, column='H2O_13.94', value=cr1000x_df_full['H2O_3'].values)
    df.insert(loc=22, column='H2Osig_13.94', value=cr1000x_df_full['H2O_sig_strgth_3'].values)
    df.insert(loc=23, column='CO2_13.94', value=cr1000x_df_full['CO2_3'].values)
    df.insert(loc=24, column='CO2sig_13.94', value=cr1000x_df_full['CO2_sig_strgth_3'].values)
    df.insert(loc=25, column='gas_diag_13.94', value=cr1000x_df_full['diag_irga_3'].values)

    df.insert(loc=26, column='Ux_6.35', value=cr1000x_df_full['Ux_1'].values)
    df.insert(loc=27, column='Uy_6.35', value=cr1000x_df_full['Uy_1'].values)
    df.insert(loc=28, column='Uz_6.35', value=cr1000x_df_full['Uz_1'].values)
    df.insert(loc=29, column='T_Sonic_6.35', value=cr1000x_df_full['Ts_1'].values)
    df.insert(loc=30, column='diagnostic_6.35', value=cr1000x_df_full['diag_sonic_1'].values)
    df.insert(loc=31, column='Pressure_6.35', value=cr1000x_df_full['cell_press_1'].values)
    df.insert(loc=32, column='H2O_6.35', value=cr1000x_df_full['H2O_1'].values)
    df.insert(loc=33, column='H2Osig_6.35', value=cr1000x_df_full['H2O_sig_strgth_1'].values)
    df.insert(loc=34, column='CO2_6.35', value=cr1000x_df_full['CO2_1'].values)
    df.insert(loc=35, column='CO2sig_6.35', value=cr1000x_df_full['CO2_sig_strgth_1'].values)
    df.insert(loc=36, column='gas_diag_6.35', value=cr1000x_df_full['diag_irga_1'].values)

    df.insert(loc=37, column='Ux_4.42', value=cr1000x_df_full['Ux_4'].values)
    df.insert(loc=38, column='Uy_4.42', value=cr1000x_df_full['Uy_4'].values)
    df.insert(loc=39, column='Uz_4.42', value=cr1000x_df_full['Uz_4'].values)
    df.insert(loc=40, column='T_Sonic_4.42', value=cr1000x_df_full['Ts_4'].values)
    df.insert(loc=41, column='diagnostic_4.42', value=cr1000x_df_full['diag_sonic_4'].values)
    df.insert(loc=42, column='Pressure_4.42', value=cr1000x_df_full['cell_press_4'].values)
    df.insert(loc=43, column='H2O_4.42', value=cr1000x_df_full['H2O_4'].values)
    df.insert(loc=44, column='H2Osig_4.42', value=cr1000x_df_full['H2O_sig_strgth_4'].values)
    df.insert(loc=45, column='CO2_4.42', value=cr1000x_df_full['CO2_4'].values)
    df.insert(loc=46, column='CO2sig_4.42', value=cr1000x_df_full['CO2_sig_strgth_4'].values)
    df.insert(loc=47, column='gas_diag_4.42', value=cr1000x_df_full['diag_irga_4'].values)

    df.apply(pd.to_numeric)
    datestr1 = f'{timestrings[0].year}{timestrings[0].month:02d}{timestrings[0].day:02d}'
    endtime_48h = timestrings[0] + pd.Timedelta(hours=48)
    datestr2 = f'{endtime_48h.year}{endtime_48h.month:02d}{endtime_48h.day:02d}'
    df.to_csv(
        FM_processed_dir + f'FMDOL_{HZ[0]}Hz_{datestr1}000000_{datestr2}000000.txt',
        header=False, index=None, sep=',')

    # # Update last_row_previous for next loop (using the **last row before reindex**)
    # last_row_previous = cr1000x_df.iloc[-1]
    del cr1000x_df, df