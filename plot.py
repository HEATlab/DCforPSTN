import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
import json

##
# \file plot.py
# \brief Methods to plot scatter plots shown in empirical sections



def plot_lsc_vs_ml(risks, lsc, ml):
    trace_lsc = go.Scatter(
        x=risks,
        y=lsc,
        name='<b>LSC-LP</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(255, 182, 193, .9)',
        ))

    trace_ml = go.Scatter(
        x=risks,
        y=ml,
        name='<b>Min-Loss DC</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(152, 0, 0, .8)',
        ))

    layout = go.Layout(
        title='<b>Dispatch Success Across Risk Levels</b>',
        titlefont=dict(
            size=24,
            color='black',
        ),
        legend=dict(
            x=0.1,
            y=1,
            # orientation="h",
            font=dict(size=20, ),
        ),
        xaxis=dict(
            title='Risk Level (α)',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )),
        yaxis=dict(
            title='Dispatch Success Rate',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )))

    data = [trace_lsc, trace_ml]
    fig = go.Figure(data=data, layout=layout)
    # py.iplot(fig, filename='test')
    py.image.save_as(fig, filename='risks_aij.png', scale=3)


##
#  \fn plot_strong_1
#  \brief plot the Actual/optimal DSC vs. LP-predicted DSC
def plot_strong_1():
    with open('result/result_compare.json', 'r') as f:
        result_compare = json.loads(f.read())

    dynamic = []
    uncertain = []

    for fname in list(result_compare.keys()):
        if fname[:7] == 'dynamic':
            dynamic.append(result_compare[fname])
        else:
            uncertain.append(result_compare[fname])

    trace0 = go.Scatter(
        x=[d[0] for d in dynamic],
        y=[d[1] for d in dynamic],
        name='<b>DC</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(255, 182, 193, .9)',
        ))

    trace1 = go.Scatter(
        x=[d[0] for d in uncertain],
        y=[d[1] for d in uncertain],
        name='<b>Not DC</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(152, 0, 0, .8)',
        ))

    layout = go.Layout(
        title='<b>Accuracy of DSC-LP Approximation</b>',
        titlefont=dict(
            size=24,
            color='black',
        ),
        legend=dict(
            x=0.1,
            y=1,
            # orientation="h",
            font=dict(size=20, ),
        ),
        xaxis=dict(
            title='LP-Predicted Degree of Strong Controllability',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )),
        yaxis=dict(
            title='Optimal Degree of Strong Controllability',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )))

    data = [trace0, trace1]
    fig = go.Figure(data=data, layout=layout)
    # py.iplot(fig, filename='test')
    py.image.save_as(fig, filename='accuracy.png', scale=3)


##
#  \fn plot_strong_2
#  \brief plot the empirical success rate vs. the LP-predicted DSC
def plot_strong_2():
    with open('result/result_success.json', 'r') as f:
        result_success = json.loads(f.read())

    dynamic = []
    uncertain = []

    for fname in list(result_success.keys()):
        if fname[:7] == 'dynamic':
            dynamic.append(result_success[fname])
        else:
            uncertain.append(result_success[fname])

    trace0 = go.Scatter(
        x=[d[0] for d in dynamic],
        y=[d[1] for d in dynamic],
        name='<b>DC</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(255, 182, 193, .9)',
        ))

    trace1 = go.Scatter(
        x=[d[0] for d in uncertain],
        y=[d[1] for d in uncertain],
        name='<b>Not DC</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(152, 0, 0, .8)',
        ))

    layout = go.Layout(
        title='<b>Fixed Decision Success Rate</b>',
        titlefont=dict(
            size=24,
            color='black',
        ),
        legend=dict(
            x=0.1,
            y=1,
            font=dict(size=20, ),
        ),
        xaxis=dict(
            title='LP-Predicted Likelihood of Strong Controllability',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )),
        yaxis=dict(
            title='Empirical Success Rate',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )))

    data = [trace0, trace1]
    fig = go.Figure(data=data, layout=layout)
    py.image.save_as(fig, filename='success.png', scale=3)


##
#  \fn plot_dynamic
#  \brief plot the empirical success rate vs. approximated DDC
def plot_dynamic():
    with open('result/dynamic_result.json', 'r') as f:
        result_dynamic = json.loads(f.read())
    with open('result/relax_result.json', 'r') as f:
        result_relax = json.loads(f.read())

    once = []
    more = []

    for fname in list(result_dynamic.keys()):
        if result_relax[fname] == 1:
            once.append(result_dynamic[fname])
        else:
            more.append(result_dynamic[fname])

    trace0 = go.Scatter(
        x=[d[0] for d in once],
        y=[d[1] for d in once],
        name='<b>Relax Once</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(255, 182, 193, .9)',
        ))

    trace1 = go.Scatter(
        x=[d[0] for d in more],
        y=[d[1] for d in more],
        name='<b>Relax Multiple Times</b>  ',
        mode='markers',
        marker=dict(
            size=13,
            color='rgba(152, 0, 0, .8)',
        ))

    layout = go.Layout(
        title='<b>Dispatch Success Rate </b>',
        titlefont=dict(
            size=24,
            color='black',
        ),
        legend=dict(
            x=0.1,
            y=1,
            font=dict(size=15, ),
        ),
        xaxis=dict(
            title='Predicted Likelihood of Dynamic Controllability',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )),
        yaxis=dict(
            title='Empirical Success Rate',
            titlefont=dict(
                size=18,
                color='black',
            ),
            zerolinewidth=2,
            tickfont=dict(
                size=15,
                color='black',
            )))

    data = [trace0, trace1]
    fig = go.Figure(data=data, layout=layout)
    py.image.save_as(fig, filename='dynamic.png', scale=3)
