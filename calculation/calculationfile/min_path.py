import os
import re

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pyvis.network import Network

from .output import Output


def read_min_path(path):
    min_path = MinPath()
    min_path.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    min_path.nsteps = int(re.findall(r"nsteps:(.+)", lines[0])[0].strip())

    mdgraph = nx.MultiDiGraph()
    nodes = [int(node) for node in lines[1].split("->")]
    mdgraph.add_nodes_from(nodes)
    
    for line in lines[2:]:
        edge = re.findall(r"(.+)->(.+):(.+)", line)

        if edge:
            node0 = int(edge[0][0].strip())
            node1 = int(edge[0][1].strip())
            eattrs = re.findall(r"PT.+?\)|CL", edge[0][2])

            for eattr in eattrs:

                if "PT" in eattr:
                    string = re.findall(r"PT(.+)\((.+)\)", eattr)
                    pt = Output()
                    pt.name = int(string[0][0].strip())
                    pt.E = float(string[0][1].strip())
                    mdgraph.add_edge(node0, node1, pt=pt)

                elif "CL" in eattr:
                    mdgraph.add_edge(node0, node1)

                else:
                    pass

    min_path.mdgraph = mdgraph

    return min_path


class MinPath:

    def __init__(self, name=None, nsteps=0, mdgraph=None):
        self.name = name
        self.nsteps = nsteps
        self.mdgraph = mdgraph

    def show(self, rpn, path="default"):
        group_eqs = {}

        for node in self.mdgraph.nodes(data=False):
            group = rpn.mdgraph.nodes[node]["eq"].group

            if group in group_eqs.keys():
                group_eqs[group].append(node)

            else:
                group_eqs[group] = [node]

        group_num = {group:len(eqs) for group, eqs in group_eqs.items()}
        group_r = {group:70*num/(2*np.pi) for group, num in group_num.items()}
        group_r_cumsum = {
            group:sum([v for k, v in group_r.items() if k <= group])
            for group in group_r.keys()
        }
        R = 3.0 * max(group_r_cumsum.values()) / np.pi

        G = nx.MultiDiGraph()
        
        for node in self.mdgraph.nodes(data=False):
            eq = rpn.mdgraph.nodes[node]["eq"]
            cmap = plt.cm.jet(eq.E/200)
            cmap255 = (cmap[0]*255, cmap[1]*255, cmap[2]*255, cmap[3])
            group = eq.group

            if group_num[group] == 1:
                r = 0
            else:
                r = group_r[group]

            T = 2*np.pi*(
                3.0*(2*group_r_cumsum[group] - group_r[group])
                /(2*np.pi*R) - 1/4
            )
            t = 2*np.pi*(
                sorted(group_eqs[group]).index(eq.name)
                /group_num[group] - 1/4
            )

            G.add_node(
                node,
                label=str(node),
                shape="dot",
                color=f"rgb{cmap255}",
                physics=False,
                title=f"G{group}",
                x=R*np.cos(T)+r*np.cos(t),
                y=R*np.sin(T)+r*np.sin(t)
            )

        for node0, node1, _, eattr in self.mdgraph.edges(data=True, keys=True):

            if eattr:
                G.add_edge(
                    node0,
                    node1,
                    physics=True,
                    title=f"PT{eattr['pt'].name}"
                )
            else:
                G.add_edge(
                    node0,
                    node1,
                    physics=True,
                    title="CL"
                )

        net = Network(height=800, width=1400, directed=True)
        net.from_nx(G)

        if path == "default":
            path = f"{self.name}.html"

        net.show(path)

        # Remove the loading bar from the html file.
        with open(path, "r") as f:
            lines = f.readlines()

        index = [
            i for i, line in enumerate(lines)
            if '<div id="loadingBar">' in line
        ]

        if index:
            lines[index[0]:index[0]+9] = []

        with open(path, "w") as f:
            f.writelines(lines)
