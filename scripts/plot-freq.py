import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# Set global default fontsizes
mpl.rcParams["legend.title_fontsize"] = 12
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["xtick.labelsize"] = 10
mpl.rcParams["ytick.labelsize"] = 10
mpl.rcParams["legend.fontsize"] = 10


def plot_freq(df_file, raw_file, color_file, output_plot):
    df = pd.read_csv(df_file, sep="\t", parse_dates=["date"])
    raw = pd.read_csv(raw_file, sep="\t", parse_dates=["date"])

    # Load color map
    colors_by_n = {}
    with open(color_file, "r", encoding="utf-8") as fh:
        for n, line in enumerate(fh):
            colors_by_n[n + 1] = line.rstrip().split("\t")

    # Add color codes by variant
    variants = sorted(df["variant"].drop_duplicates().values)
    color_by_variant = dict(zip(variants, colors_by_n[len(variants)]))
    df["variant_color"] = df["variant"].map(color_by_variant)

    fig = sns.FacetGrid(
        data=df,
        col="location",
        hue="variant",
        col_wrap=4,
        palette=color_by_variant,
        hue_order=list(color_by_variant.keys())
    ).set_axis_labels(
        x_var="Date",
        y_var="Frequency"
    ).set_titles(
        col_template="{col_name}"
    ).tick_params(
        axis="x",
        rotation=45
    )

    # Plot the CIs
    def plot_with_ci(data, **kwargs):
        sns.lineplot(data=data, x="date", y="median", **kwargs)
        plt.fill_between(data["date"], data["HDI_95_lower"], data["HDI_95_upper"], alpha=0.2, color=kwargs["color"])
    fig.map_dataframe(plot_with_ci)

    # Plot the weekly_raw_freq
    def plot_with_raw(data, **kwargs):
        raw_data = raw[(raw["location"] == data["location"].iloc[0]) & (raw["variant"] == data["variant"].iloc[0])]
        sns.scatterplot(data=raw_data, x="date", y="raw_freq", **kwargs, s=10, alpha = 0.7, legend=False)
    fig.map_dataframe(plot_with_raw)

    fig.add_legend(title="Variants")
    sns.move_legend(fig, loc="center right")
    fig.savefig(output_plot, bbox_inches="tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Freq plots by location and variant.")
    parser.add_argument("-i", "--input_freq", type=str, required=True, help="Parsed MLR site freq TSV file (<model>_freq.tsv)")
    parser.add_argument("-r", "--input_raw", type=str, required=True, help="Path to weekly raw sequence TSV file")
    parser.add_argument("-c", "--colors", type=str, required=True, help="Path to Nextstrain color scheme (configs/color_schemes.tsv)]")
    parser.add_argument("-o", "--output", type=str, required=True, help="Site frequency by location plot PDF")
    args = parser.parse_args()
    plot_freq(args.input_freq, args.input_raw, args.colors, args.output)