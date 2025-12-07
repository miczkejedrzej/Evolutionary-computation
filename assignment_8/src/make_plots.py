import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === CONFIG ===
directory = "../results/"
separator = ";"                         # CSV separator
xcol, ycol = "x", "y"                   # Column names


# === Load all CSV files ===
for dataset in ['A', 'B']:
    files = [f for f in os.listdir(directory) if f.endswith(".csv") and f.startswith(dataset)]

    if not files:
        print("No CSV files found in directory:", directory)
        exit()

    datasets = {}
    all_x, all_y = [], []

    for filename in files:
        path = os.path.join(directory, filename)

        df = pd.read_csv(path, sep=separator, header=None, names=[xcol, ycol])
        datasets[filename] = df

        all_x.extend(df[xcol].tolist())
        all_y.extend(df[ycol].tolist())


    # === Establish common axis ranges ===
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)


    # === Plot each file in its own figure ===
    for filename, df in datasets.items():
        name = os.path.splitext(filename)[0]  # remove ".csv"

        plt.figure(figsize=(8, 5))

        # ----- Scatter plot -----
        plt.scatter(df[xcol], df[ycol], alpha=0.3, s=10, label="Data")

        # ----- Linear regression -----
        x = df[xcol].values
        y = df[ycol].values

        # Compute slope m and intercept b for y = m*x + b
        m, b = np.polyfit(x, y, 1)

        # Generate line endpoints over global x-range
        line_x = np.array([xmin, xmax])
        line_y = m * line_x + b

        # ----- Correlation coefficient -----
        r = np.corrcoef(x, y)[0, 1]   # Pearson's r

        # Plot regression line
        plt.plot(
            line_x, line_y,
            color="red",
            linewidth=2,
            label=f"Linear Regression\nCorrelation Coefficient = {r:.4f}"
        )

        # ----- Titles and labels -----
        plt.title(name)
        plt.xlabel("Objective value")
        plt.ylabel("(Average) Similarity [%]")
        plt.grid(True)
        plt.legend()

        # ----- Shared axis limits -----
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

        plt.tight_layout()
        
        # Save image
        plt.savefig(f"../visualizations/{name}.png")

plt.show()
