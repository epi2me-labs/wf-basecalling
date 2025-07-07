#!/usr/bin/env python
"""Create tables for the report."""
from bokeh.models import Title
from ezcharts.plots.distribution import histplot
import pandas as pd

# PLOTS

# The SeqSummary from ezcharts.components.fastcat cannot be used.
# It groups data into bins, but from the real time analysis output
# the input data is already grouped into bins.
# Use weights of histplot for y axis.


def read_quality_plot(seq_summary, min_qual=4, max_qual=30, title='Read quality'):
    """Create read quality summary plot."""
    df = pd.DataFrame.from_dict(seq_summary['qual'].items())
    df.columns = ['mean_quality', 'counts']
    df['mean_quality'] = df['mean_quality'].astype('float')
    plt = histplot(
        data=df['mean_quality'],
        bins=len(df),
        weights=list(df['counts'])
        )
    plt._fig.add_layout(
        Title(text=title, text_font_size="1.5em"),
        'above'
    )
    plt._fig.xaxis.axis_label = "Quality score"
    plt._fig.yaxis.axis_label = "Number of reads"
    plt._fig.x_range.start = min_qual
    plt._fig.x_range.end = max_qual
    return plt


def read_length_plot(seq_summary, title='Read length'):
    """Create a read length plot."""
    df = pd.DataFrame.from_dict(seq_summary['len'].items())
    df.columns = ['read_length', 'counts']
    df['read_length'] = df['read_length'].astype('uint64')
    df['read_length'] = df['read_length'] / 1000
    plt = histplot(
        data=df['read_length'],
        bins=len(df),
        weights=list(df['counts']))
    plt._fig.add_layout(
        Title(text=title, text_font_size="1.5em"),
        'above'
    )
    plt._fig.x_range.start = 0
    plt._fig.xaxis.axis_label = "Read length / kb"
    plt._fig.yaxis.axis_label = "Number of reads"
    return plt
