"""Make a HTML file for the docs, plus copy all the files."""


import os
import shutil
import sys


sys.stderr = sys.stdout = open(snakemake.log[0], "w")


docsdir = snakemake.output.docsdir
charts = snakemake.params.charts
all_inputs = snakemake.input


def extract_final_values(d):
    """Get list of final values in nested dict."""
    values = []
    for key, value in d.items():
        if isinstance(value, dict):
            values.extend(extract_final_values(value))
        else:
            values.append(value)
    return values


if not os.path.isdir(docsdir):
    os.makedirs(docsdir, exist_ok=True)

chartbases = set([])
assert set(extract_final_values(charts)).issubset(all_inputs)
for chart in all_inputs:
    chartbase = os.path.basename(chart)
    if chartbase in chartbases:
        raise ValueError(f"Duplicate {chartbase=} {chart=}")
    chartbases.add(chartbase)
    shutil.copy(chart, os.path.join(docsdir, chartbase))


def dict_to_html(d, level=1):
    html_content = ""
    for key, value in d.items():
        if isinstance(value, dict):
            html_content += f"<h{level}>{key}</h{level}>\n"
            html_content += dict_to_html(value, level + 1)
        else:
            value_base = os.path.basename(value)
            html_content += f'<div><a href="{value_base}">{key}</a></div>'
    return html_content


title = snakemake.params.title
heading = snakemake.params.heading
html_document = (
    f"<!DOCTYPE html>\n<html>\n<head>\n<title>{title}</title>\n</head>\n<body>\n"
    f"<h1>{title}</h1>\n<div>{heading}</div>\n"
)
html_document += dict_to_html(charts, level=2)
html_document += "</body>\n</html>"
    
with open(os.path.join(docsdir, "index.html"), "w") as f:
    f.write(html_document)
