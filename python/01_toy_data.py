# Load and inspect the toy dataset

import os
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(script_dir, "..", "data", "toy_microbiome_metabolite_outcome.csv")
df = pd.read_csv(csv_path)

# View the data
print(df)

# Summary statistics
print(df.describe())
