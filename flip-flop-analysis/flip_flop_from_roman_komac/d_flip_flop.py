import pandas as pd
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def plot_results():
    df1 = pd.read_csv('Concentrations_average_1.csv')

    fig1 = make_subplots(rows=2, cols=2)

    fig1.add_trace(go.Scatter(name='activatory proteins', x=df1['time(seconds)'], y=df1['# activatory proteins']), row=1, col=1)
    fig1.add_trace(go.Scatter(name='repressory proteins', x=df1['time(seconds)'], y=df1['# repressory proteins']), row=1, col=2)
    fig1.add_trace(go.Scatter(name='Q proteins', x=df1['time(seconds)'], y=df1['# Q proteins']), row=2, col=1)
    fig1.add_trace(go.Scatter(name='Qc proteins', x=df1['time(seconds)'], y=df1['# Qc proteins']), row=2, col=2)

    plotly.offline.plot(fig1, filename='concentrations_result_1.html')

    df2 = pd.read_csv('Concentrations_average_2.csv')

    fig2 = make_subplots(rows=1, cols=3)

    fig2.add_trace(go.Scatter(name='repressory proteins', x=df2['time(seconds)'], y=df2['# repressory proteins']), row=1, col=1)
    fig2.add_trace(go.Scatter(name='Q proteins', x=df2['time(seconds)'], y=df2['# Q proteins']), row=1, col=2)
    fig2.add_trace(go.Scatter(name='Qc proteins', x=df2['time(seconds)'], y=df2['# Qc proteins']), row=1, col=3)

    plotly.offline.plot(fig2, filename='concentrations_result_2.html')


if __name__ == '__main__':
    plot_results()
