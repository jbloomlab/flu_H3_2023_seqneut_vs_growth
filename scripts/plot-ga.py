from matplotlib import pyplot as plt
import pandas as pd
import argparse
import math
            
def plot_ga(input_file, location, virus, color_file, out_var, out_loc):
    # Read in ga data and remove location "hierarchical"
    df = pd.read_csv(input_file, sep="\t").query("location != 'hierarchical'")

    # Load color map
    colors_by_n = {}
    with open(color_file, "r", encoding="utf-8") as fh:
        for n, line in enumerate(fh):
            colors_by_n[n + 1] = line.rstrip().split("\t")

    # Add color codes by location
    locations = sorted(df["location"].drop_duplicates().values)
    color_by_location = dict(zip(locations, colors_by_n[len(locations)]))
    df["location_color"] = df["location"].map(color_by_location)
    
    # Add color codes by variant
    variants = sorted(df["variant"].drop_duplicates().values)
    color_by_variant = dict(zip(variants, colors_by_n[len(variants)]))
    df["variant_color"] = df["variant"].map(color_by_variant)

    # Plot GA by location
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=300)

    for location_i, location_df in df.groupby("location"):
        color = color_by_location[location_i]
    
        ax.plot(
            location_df["median"],
            location_df["variant"],
            "o",
            alpha=0.5,
            color=color,
            label=location_i
            )

    ax.axvline(
        x=1.0,
        color="black",
        linestyle="dashed",
        zorder=-10
        )

    ax.legend(
        title=location.capitalize(),
        frameon=False,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        borderaxespad=0,
        title_fontsize="large",
        labelspacing=0.6,
        handleheight=0.6
        )

    ax.set_title(f"{virus} Growth Advantage (GA) by {location.capitalize()}", fontsize=12)
    ax.set_xlabel("Median GA", fontsize=12)
    ax.set_ylabel("Subclade", fontsize=12)

    fig.savefig(out_loc, bbox_inches="tight")
    
    # Plot GA by subclade
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=300)

    for variant_i, variant_df in df.groupby("variant"):
        color = color_by_variant[variant_i]
    
        ax.plot(
            variant_df["median"],
            variant_df["location"],
            "o",
            alpha=0.5,
            color=color,
            label=variant_i
            )

    ax.axvline(
        x=1.0,
        color="black",
        linestyle="dashed",
        zorder=-10
        )

    ax.legend(
        title="Subclade",
        frameon=False,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0,
        ncol=2,
        title_fontsize="large",
        labelspacing=1,
        handleheight=1
        )

    ax.set_title(f"{virus} Growth Advantage (GA) by Subclade", fontsize=12)
    ax.set_xlabel("Median GA", fontsize=12)
    ax.set_ylabel(location.capitalize(), fontsize=12)

    fig.savefig(out_var, bbox_inches="tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot GA plots by location and variant.")
    parser.add_argument("-i", "--input_ga", type=str, help="Parsed MLR growth advantage file (<name>_ga.tsv)")
    parser.add_argument("-l", "--location", help="String specifying location ['region', 'country']")
    parser.add_argument("-v", "--virus", help="Virus type ['H3N2', 'H1N1pdm', 'B_Vic']")
    parser.add_argument("-c", "--colors", help="Path to Nextstrain color scheme (color_schemes.tsv)]")
    parser.add_argument("-ov", "--out_variant", help="GA by variant PDF")
    parser.add_argument("-ol", "--out_location", help="GA by location PDF")
    args = parser.parse_args()
    plot_ga(args.input_ga, args.location, args.virus, args.colors, args.out_variant, args.out_location)
