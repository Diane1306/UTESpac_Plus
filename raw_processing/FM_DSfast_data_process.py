import os
import pandas as pd

c_path = os.getcwd()
if 'diane_wt' in c_path:
    box_path = '/Users/diane_wt/Library/CloudStorage/Box-Box'
elif 'admin-dianew68' in c_path:
    box_path = 'C:/Users/admin-dianew68/Box'

FM_dir = box_path + '/Lab Library/French Meadows/Summer2025/data/20250527_ds/cr3000/ascii/' # FM DS Tower data dir
FM_processed_dir = box_path + "/Lab Library/French Meadows/Summer2025/data/processed/FM_DS_10Hz/" # processed data dir

startid = 1380 # **** change this index corresponds to a start date
endid = 4787   # **** change this index corresponds to an end date
HZ = [10]

last_row_previous = None  # initialize outside loop
for fi in range(startid, endid, 48):
    cr3000_files = [f"{FM_dir}TOA5_5500.FastResponse_{fi + ii}.dat" for ii in range(48)]
    existing_files = [f for f in cr3000_files if os.path.exists(f)]

    if not existing_files:
        print(f"⚠️ Skipping batch starting at file index {fi}: no files found.")
        continue
    try:
        cr3000_df = pd.concat([
            pd.read_csv(f, skiprows=[0, 2, 3], index_col=[0], na_values=['NaN', 'NAN'], parse_dates=True)
            for f in existing_files
        ], axis=0)
    except Exception as e:
        print(f"❌ Failed to load batch starting at {fi} due to error: {e}")
        continue

    cr3000_df = cr3000_df[~cr3000_df.index.duplicated(keep='first')]
    timestrings = pd.to_datetime(cr3000_df.index, format='mixed')

    # Create full time index covering 48 hours at 0.1 s step
    start_time = pd.to_datetime(timestrings[0].floor('D'))  # midnight of first day
    end_time = start_time + pd.Timedelta(hours=48)
    full_time_index = pd.date_range(start=start_time, end=end_time, freq='100ms', inclusive='left')
    # Reindex your dataframe to this full time index
    cr3000_df_full = cr3000_df.reindex(full_time_index)
    # Continue your processing pipeline using cr3000_df_full instead of cr3000_df

    # Only from second loop onwards: fill first row with last_row_previous if missing
    if last_row_previous is not None:
        if cr3000_df_full.iloc[0].isna().all():
            cr3000_df_full.iloc[0] = last_row_previous.values

    timestrings = cr3000_df_full.index  # update timestrings to full time index

    print(f"Expected rows: {len(full_time_index)}, Actual rows before reindex: {len(cr3000_df)}, and after reindex: {len(cr3000_df_full)}")

    df = pd.DataFrame({})
    df.insert(loc=0, column='year', value=[int(stamp.year) for stamp in timestrings])
    df.insert(loc=1, column='day', value=[int(stamp.day_of_year) for stamp in timestrings])
    df.insert(loc=2, column='HM', value=[int(f"{stamp.hour:02d}{stamp.minute:02d}") for stamp in timestrings])
    df.insert(loc=3, column='second', value=[f"{stamp.second + stamp.microsecond/1000000:.1f}" for stamp in timestrings]) # *****Change here if data is in 20 Hz*****

    df.insert(loc=4, column='Ux_6.83', value=cr3000_df_full['Ux_2_0'].values)
    df.insert(loc=5, column='Uy_6.83', value=cr3000_df_full['Uy_2_0'].values)
    df.insert(loc=6, column='Uz_6.83', value=cr3000_df_full['Uz_2_0'].values)
    df.insert(loc=7, column='T_Sonic_6.83', value=cr3000_df_full['Ts_2_0'].values)
    df.insert(loc=8, column='diagnostic_6.83', value=cr3000_df_full['diag_sonic_2_0'].values)

    df.insert(loc=9, column='Ux_4.13', value=cr3000_df_full['Ux_1_0'].values)
    df.insert(loc=10, column='Uy_4.13', value=cr3000_df_full['Uy_1_0'].values)
    df.insert(loc=11, column='Uz_4.13', value=cr3000_df_full['Uz_1_0'].values)
    df.insert(loc=12, column='T_Sonic_4.13', value=cr3000_df_full['Ts_1_0'].values)
    df.insert(loc=13, column='diagnostic_4.13', value=cr3000_df_full['diag_sonic_1_0'].values)

    df.apply(pd.to_numeric)
    datestr1 = f'{timestrings[0].year}{timestrings[0].month:02d}{timestrings[0].day:02d}'
    endtime_48h = timestrings[0] + pd.Timedelta(hours=48)
    datestr2 = f'{endtime_48h.year}{endtime_48h.month:02d}{endtime_48h.day:02d}'
    df.to_csv(
        FM_processed_dir + f'FMDS_{HZ[0]}Hz_{datestr1}000000_{datestr2}000000.txt',
        header=False, index=None, sep=',')

    # Update last_row_previous for next loop (using the **last row before reindex**)
    last_row_previous = cr3000_df.iloc[-1]
    del cr3000_df, df