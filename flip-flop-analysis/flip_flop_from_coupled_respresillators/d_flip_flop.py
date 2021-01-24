import pandas as pd
import plotly
import plotly.express as px


def plot_results():
    df_ai = pd.read_csv('AI_internal_ALL.csv')
    fig_ai = px.line(df_ai, x='time(seconds)', y='internal_AI', title='Plot')
    plotly.offline.plot(fig_ai, filename='ai_internal_result.html')

    df_q = pd.read_csv('chemicalFields.csv')
    fig_q = px.line(df_q, x='time(seconds)', y='qFieldAvg', title='Plot')
    plotly.offline.plot(fig_q, filename='q_chem_field_result.html')

    df_qc = pd.read_csv('chemicalFields.csv')
    fig_qc = px.line(df_qc, x='time(seconds)', y='qcFieldAvg', title='Plot')
    plotly.offline.plot(fig_qc, filename='qc_chem_field_result.html')

    df_lacI = pd.read_csv('lacI_ALL.csv')
    fig_lacI = px.line(df_lacI, x='time(seconds)', y='lacI_mRNA', title='Plot')
    plotly.offline.plot(fig_lacI, filename='lacI_result.html')


if __name__ == '__main__':
    plot_results()
