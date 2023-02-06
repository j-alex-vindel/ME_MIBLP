import sys
import pandas as pd
import ast
from typing import Type,List
import matplotlib.pyplot as plt
import seaborn as sns


bacteria = Type[str]
graphs = List[str]

def make_plot(name:bacteria=None,choice:graphs=None):

    s1 = pd.read_csv(f"../Results/EQ/Methods_K1_{name}.csv")
    s2 = pd.read_csv(f"../Results/EQ/Methods_K2_{name}.csv")
    s3 = pd.read_csv(f"../Results/EQ/Methods_K3_{name}.csv")

    st = pd.concat([s1,s2,s3],sort=False).round(decimals=5)

    st['K'] = st.apply(lambda row: len(ast.literal_eval(row.Strategy)),axis=1)

    if 't' in choice:
        timeplot = sns.catplot(x='Method',y='Time',hue='K',kind='point',data=st,palette='flare',linestyles=['dashed','dotted','dashdot'],markers=['^','o','s'])
        plt.savefig(f"../Results/Graphs/Time_lines_{name}.png")
        sns.catplot(data=st,x='Method',y='Time',hue='K',kind='bar',palette='flare')
        plt.savefig(f"../Results/Graphs/Time_bar_{name}.png")

    if 'b' in choice:
        bp = sns.catplot(data=st, x="Biomass", y="Method", hue='K',kind='swarm',palette='flare')
        plt.savefig(f"../Results/Graphs/Biomass_{name}.png")

    if 'c' in choice:
        cp = sns.catplot(data=st, x="Chemical", y="Method", hue='K',kind='swarm',palette='flare')
        plt.savefig(f"../Results/Graphs/Chemical_{name}.png")

    if 'tm' in choice:
        sns.catplot(data=st,x='K',y='Time',hue='Method',kind='bar',palette='flare')
        plt.savefig(f"../Results/Graphs/Time_bymethod_{name}.png")

if __name__ =='__main__':
    bacterias = ['iJO1366','iAF1260','iJR904']
    choice = 'tm'
    for bact in bacterias:
        make_plot(name=bact,choice=choice)