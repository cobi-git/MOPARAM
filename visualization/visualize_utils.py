def change_column_name(df):
    df.columns = [a.replace('_test', '').replace('omics_combination_', '').replace('_', ' ').capitalize() for a in df.columns]