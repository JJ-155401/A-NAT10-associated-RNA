import pandas as pd

sc = pd.read_csv('results/single_cell_analysis.csv')
print('=== Correct cell counts from single_cell_analysis.csv ===')
total = int(sc['n_ctrl'].sum() + sc['n_ad'].sum())
print(f'Total filtered cells: {total}')
for _, row in sc.iterrows():
    n = int(row['n_ctrl'] + row['n_ad'])
    pct = 100 * n / total
    print(f"  {row['cell_type']}: {n} ({pct:.1f}%)")
