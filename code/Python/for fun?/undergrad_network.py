import pandas as pd
import igraph

df = pd.read_excel("../../../data/my_undergrad_data/2017-18_bugs.xlsx", sheet_name="Master", header=0,)

df = df[["Stream", "Reach", "Treatment", "Meter", "Count", "Order", "Family", "Genus"]]

df = df.groupby(by=["Stream", "Reach", "Treatment", "Family", "Genus"], as_index=False).sum()

df = df.pivot(columns="Family", values="Count")

df = df.fillna(value=0)

df = df.astype(bool).astype(int) #turn all values into either 0 or 1

cooccurance_matrix = df.T.dot(df)

print(cooccurance_matrix.head(6))


