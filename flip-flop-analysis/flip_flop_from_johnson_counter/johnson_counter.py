import pandas as pd
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def plot_results():
    # first flip flop outputs
    df1 = pd.read_csv('flip_flop1_concentrations_average.csv')

    fig1 = make_subplots(rows=2, cols=2)

    fig1.add_trace(go.Scatter(name='activatory proteins', x=df1['time(seconds)'], y=df1['activatory proteins(h)']), row=1, col=1)
    fig1.add_trace(go.Scatter(name='repressory proteins', x=df1['time(seconds)'], y=df1['repressory proteins(i)']), row=1, col=2)
    fig1.add_trace(go.Scatter(name='Q proteins', x=df1['time(seconds)'], y=df1['Q proteins(q)']), row=2, col=1)
    fig1.add_trace(go.Scatter(name='Qc proteins', x=df1['time(seconds)'], y=df1['Qc proteins(qc)']), row=2, col=2)

    plotly.offline.plot(fig1, filename='concentrations_flip_flop1.html')

    # second flip flop outputs
    df2 = pd.read_csv('flip_flop2_concentrations_average.csv')

    fig2 = make_subplots(rows=2, cols=2)

    fig2.add_trace(go.Scatter(name='activatory proteins', x=df2['time(seconds)'], y=df2['activatory proteins(h)']), row=1, col=1)
    fig2.add_trace(go.Scatter(name='repressory proteins', x=df2['time(seconds)'], y=df2['repressory proteins(i)']), row=1, col=2)
    fig2.add_trace(go.Scatter(name='Q proteins', x=df2['time(seconds)'], y=df2['Q proteins(q)']), row=2, col=1)
    fig2.add_trace(go.Scatter(name='Qc proteins', x=df2['time(seconds)'], y=df2['Qc proteins(qc)']), row=2, col=2)

    plotly.offline.plot(fig2, filename='concentrations_flip_flop2.html')

    # third flip flop outputs
    df3 = pd.read_csv('flip_flop3_concentrations_average.csv')

    fig3 = make_subplots(rows=2, cols=2)

    fig3.add_trace(go.Scatter(name='activatory proteins', x=df3['time(seconds)'], y=df3['activatory proteins(h)']), row=1, col=1)
    fig3.add_trace(go.Scatter(name='repressory proteins', x=df3['time(seconds)'], y=df3['repressory proteins(i)']), row=1, col=2)
    fig3.add_trace(go.Scatter(name='Q proteins', x=df3['time(seconds)'], y=df3['Q proteins(q)']), row=2, col=1)
    fig3.add_trace(go.Scatter(name='Qc proteins', x=df3['time(seconds)'], y=df3['Qc proteins(qc)']), row=2, col=2)

    plotly.offline.plot(fig3, filename='concentrations_flip_flop3.html')

    # third flip flop outputs (longer simulation)
    df3_longer = pd.read_csv('flip_flop3_concentrations_average_longer.csv')

    fig3_longer = make_subplots(rows=2, cols=2)

    fig3_longer.add_trace(go.Scatter(name='activatory proteins', x=df3_longer['time(seconds)'], y=df3_longer['activatory proteins(h)']), row=1, col=1)
    fig3_longer.add_trace(go.Scatter(name='repressory proteins', x=df3_longer['time(seconds)'], y=df3_longer['repressory proteins(i)']), row=1, col=2)
    fig3_longer.add_trace(go.Scatter(name='Q proteins', x=df3_longer['time(seconds)'], y=df3_longer['Q proteins(q)']), row=2, col=1)
    fig3_longer.add_trace(go.Scatter(name='Qc proteins', x=df3_longer['time(seconds)'], y=df3_longer['Qc proteins(qc)']), row=2, col=2)

    plotly.offline.plot(fig3_longer, filename='concentrations_flip_flop3_longer.html')


if __name__ == '__main__':
    plot_results()
