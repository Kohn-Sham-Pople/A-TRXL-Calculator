from io import BytesIO
import plotly.graph_objects as go
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def interactive(x,series,xlabel,title):
    fig=go.Figure()
    for label,y in series.items(): fig.add_scatter(x=x,y=y,name=label,mode="lines")
    fig.update_layout(title=title,xaxis_title=xlabel,yaxis_title="Pair weight / contrast (electrons²)",template="plotly_white")
    return fig

def png(x,series,xlabel,title):
    fig,ax=plt.subplots(figsize=(8,4.5),layout="constrained")
    for label,y in series.items(): ax.plot(x,y,label=label)
    ax.set(xlabel=xlabel,ylabel="Pair weight / contrast (electrons²)",title=title)
    ax.legend();ax.grid(alpha=.2)
    out=BytesIO();fig.savefig(out,format="png",dpi=180);plt.close(fig)
    return out.getvalue()
