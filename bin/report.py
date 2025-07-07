#!/usr/bin/env python
"""Create workflow report."""
import argparse
import json

from bokeh.models import Title
from dominate.tags import p
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
import pandas as pd
from report_utils import read_length_plot, read_quality_plot
from ezcharts.util import get_named_logger  # noqa: ABS101


THEME = 'epi2melabs'


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Basecalling report", "wf-basecalling",
        args.params, args.versions,
        args.workflow_version,)

    # Create statistics from the histogram json generated from
    #  the workflow, summarising the results from bamstats.
    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            with open(args.stats[0]) as f:
                datas = json.load(f)
                tabs = Tabs()
                total_reads = {}
                for sample_id, data in sorted(datas.items()):
                    with tabs.add_tab(sample_id):
                        total_reads[sample_id] = data['total_reads']
                        if data['total_reads'] == 0:
                            p("No reads called.")
                            continue
                        with Grid(columns=2):
                            EZChart(
                                read_quality_plot(data), THEME)
                            EZChart(read_length_plot(data), THEME)
                with tabs.add_tab('total'):
                    with Grid(columns=1):  # total read counts per sample
                        df_stats = pd.DataFrame.from_dict(total_reads.items())
                        df_stats.columns = ['Sample_name', 'Number of reads']
                        plt = ezc.barplot(
                            data=df_stats, x='Sample_name', y='Number of reads')
                        plt._fig.add_layout(
                            Title(
                                text="Number of reads per sample.",
                                text_font_size="1.5em"
                            ),
                            'above'
                        )
                        EZChart(plt, THEME)

    # If pairing rates are provided, show them.
    if args.pairings:
        with report.add_section("Pairing summary", "Pairing summary"):
            with open(args.pairings[0]) as f:
                # Load data
                data = pd.read_csv(f)
                if data.empty:
                    p("No reads called.")
                else:
                    # Make summary
                    data_sum = data\
                        .drop(columns=['Filename'])\
                        .sum()\
                        .to_frame()\
                        .T
                    data_sum['Pairing rate'] = data_sum['Paired'] / data_sum['Simplex']
                    data_sum['Pairing rate'] = data_sum['Pairing rate'].round(4)
                    DataTable.from_pandas(
                        data_sum, use_index=False, export=True,
                        file_name=(
                            f'{args.sample_name}-wf-basecalling-duplex-summary'))
            p(
                'Simplex: the number of initial reads.',
                'Paired: the number of simplex reads belonging to a pair.',
                'Duplex: the number of duplex reads.',
            )

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--pairings", nargs='*', help="Pairing per-chunk stats.", required=False)
    parser.add_argument(
        "--sample_name", required=True, help="Sample name.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--workflow_version", default='unknown',
        help="Workflow version")
    return parser


if __name__ == "__main__":
    parser = argparser()
    main(parser.parse_args())
