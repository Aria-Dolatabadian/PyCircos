#https://github.com/ponnhide/pyCircos

#install python-circos
#Tutorial 1
import matplotlib.pyplot as plt
from pycircos import Gcircle
from pycircos import Garc
circle = Gcircle()
with open("example_data_chromosome_general.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",")
        name   = line[0]
        length = int(line[-1])
        arc    = Garc(arc_id=name, size=length, interspace=3, raxis_range=(950,1000), labelposition=60, label_visible=True)
        circle.add_garc(arc)

circle.set_garcs()

#cytoband
import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27",
               "acen":"#D82322"}

arcdata_dict = collections.defaultdict(dict)
with open("example_data_chromosome_cytoband.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]
        start = int(line[1])-1
        width = int(line[2])-(int(line[1])-1)
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = []
            arcdata_dict[name]["colors"]    = []
        arcdata_dict[name]["positions"].append(start)
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["colors"].append(color_dict[line[-1]])

for key in arcdata_dict:
    circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"],
                   width=arcdata_dict[key]["widths"], raxis_range=[950,1000], facecolor=arcdata_dict[key]["colors"])

# scatter plot
values_all = []
arcdata_dict = collections.defaultdict(dict)
with open("example_data_point.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        name = line[0]
        start = int(line[1]) - 1
        end = int(line[2])
        mid = (start + end) / 2
        value = float(line[-1])
        values_all.append(value)
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["values"] = []
        arcdata_dict[name]["positions"].append(mid)
        arcdata_dict[name]["values"].append(value)

vmin, vmax = min(values_all), max(values_all)
for key in arcdata_dict:
    circle.scatterplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                       rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)], raxis_range=(860, 940),
                       facecolor="orangered", spine=True)

# line plot
values_all = []
arcdata_dict = collections.defaultdict(dict)
with open("example_data_point.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        name = line[0]
        start = int(line[1]) - 1
        end = int(line[2])
        mid = (start + end) / 2
        value = float(line[-1])
        values_all.append(value)
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["values"] = []
        arcdata_dict[name]["positions"].append(mid)
        arcdata_dict[name]["values"].append(value)

vmin, vmax = min(values_all), max(values_all)
for key in arcdata_dict:
    circle.lineplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                    rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)], raxis_range=(770, 850),
                    linecolor="royalblue", spine=False)

#bar plot
values_all   = []
arcdata_dict = collections.defaultdict(dict)
with open("example_data_barplot.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]
        start = int(line[1])-1
        end   = int(line[2])
        width = end-start
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = []
            arcdata_dict[name]["values"]    = []
        arcdata_dict[name]["positions"].append(start)
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all)
for key in arcdata_dict:
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[680,760], facecolor="y", spine=True)

#heatmap
values_all   = []
arcdata_dict = collections.defaultdict(dict)
with open("example_data_rect_gradual.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]
        start = int(line[1])-1
        end   = int(line[2])
        width = end-start
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = []
            arcdata_dict[name]["values"]    = []
        arcdata_dict[name]["positions"].append(start)
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all)
for key in arcdata_dict:
    circle.heatmap(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                   width=arcdata_dict[key]["widths"], raxis_range=[630,670], vmin=vmin, vmax=vmax,
                   cmap=plt.cm.viridis)

#linkplot
#heatmap
values_all   = []
arcdata_dict = collections.defaultdict(dict)
with open("example_data_links.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name1  = line[0]
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 630)
        destination = (name2, start2, end2, 630)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)


plt.show() #or
circle.figure.savefig("Fig1.jpg")






# Toturial 2

from Bio import SeqIO
import matplotlib.pyplot as plt
from pycircos import Gcircle
from pycircos import Garc
record = SeqIO.read("NC_000913.gbk", format="genbank")
garc   = Garc(arc_id="NC_000913", record=record, interspace=0, linewidth=0,
              facecolor="#FFFFFF00", raxis_range=(0,10),
              label="Escherichia coli str. K-12 substr. MG1655", label_visible=True)

gcircle = Gcircle()
gcircle.add_garc(garc)
gcircle.set_garcs()

#Plot CDS
plus_CDS  = []
minus_CDS = []
for feat in garc.record.features:
    if feat.type == "CDS" and feat.strand >= 0:
        plus_CDS.append(feat)
    elif feat.strand == -1:
        minus_CDS.append(feat)
gcircle.featureplot("NC_000913", source=plus_CDS,  raxis_range=(700,780), facecolor="tomato")
gcircle.featureplot("NC_000913", source=minus_CDS, raxis_range=(780,860), facecolor="cornflowerblue")

#Plot GCskew
import copy
skews = garc.calc_nnskew(n1="G", n2="C")
positive_skews=copy.deepcopy(skews)
positive_skews[skews<0]=0
negative_skews=copy.deepcopy(skews)
negative_skews[skews>=0]=0
gcircle.fillplot("NC_000913", positive_skews, rlim=(min(skews),max(skews)), base_value=0, raxis_range=(400,700), facecolor="r")
gcircle.fillplot("NC_000913", negative_skews, rlim=(min(skews),max(skews)), base_value=0, raxis_range=(400,700), facecolor="b")

plt.show() #or
gcircle.figure.savefig("Fig2.jpg")






#Tutorial 3
from Bio import SeqIO
import matplotlib.pyplot as plt
from pycircos import Gcircle
from pycircos import Garc

record = SeqIO.read("NC_000913.gbk", format="genbank")
garc   = Garc(arc_id="NC_000913.3", record=record, interspace=0, linewidth=0,
              facecolor="#FFFFFF00", raxis_range=(0,10), label_visible=False)

gcircle = Gcircle()
gcircle.add_garc(garc)
gcircle.set_garcs()

#calc CDS density
plus_CDS  = []
minus_CDS = []
for feat in garc.record.features:
    if feat.type == "CDS" and feat.strand >= 0:
        plus_CDS.append((feat.location.parts[0].start, feat.location.parts[-1].end))
    elif feat.strand == -1:
        minus_CDS.append((feat.location.parts[-1].start, feat.location.parts[0].end))
plus_density  = garc.calc_density(plus_CDS, window_size=10000)
minus_density = garc.calc_density(minus_CDS, window_size=10000)
gcircle.heatmap("NC_000913.3", plus_density,  raxis_range=(700,780), cmap=plt.cm.Reds)
gcircle.heatmap("NC_000913.3", minus_density, raxis_range=(780,860), cmap=plt.cm.Blues)

# cord plot
import collections

chord_dict = collections.defaultdict(list)
with open("segdup.txt", "r") as f:
    for line in f:
        line = line.rstrip().split("\t")
        chord_dict[line[0]].append((line[1], int(line[2]), int(line[3]), 700))

for key in chord_dict:
    gcircle.chord_plot(chord_dict[key][0], chord_dict[key][1], facecolor="#ff8c0080")

plt.show() #or
gcircle.figure.savefig("Fig3.jpg")




#Tutorial 4

import pycircos
Tarc    = pycircos.Tarc
Tcircle = pycircos.Tcircle

import os
print(os.getcwd())
tarc    = Tarc(tree="kegg.nwk", format="newick", interspace=1)
tcircle = Tcircle(figsize=(12,12))

#Set colors for terminal clades
phyla_names  = ["Actinobacteria","Aquificae","Bacteroidetes", "Chlamydiae","Chlorobi","Chloroflexi","Crenarchaeota",
                "Cyanobacteria","Euryarchaeota","Firmicutes", "Spirochaetes", "Proteobacteria", "Tenericutes","Thermi","Thermotogae","Other"]
phyla_colors = ["#9ACD32", "#EE6A50", "#87CEFA", "#FFC125", "#D15FEE", "#8DEEEE", "#800000","#006400", "#800080",
                "#808080", "#B0171F", "#B0171F", "#191970", "#7B68EE", "#00CD00", "#000000"]
pylum_color_dict = dict(zip(phyla_names, phyla_colors))

cladevisual_dict  = {}
phylum_clade_dict = {}
with open("tippoint_attr.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        if line[1] in pylum_color_dict:
            cladevisual_dict[line[0]] = {"color":pylum_color_dict[line[1]], "size":12, "linewidth":0.1, "edgecolor":"#303030"}
        else:
            cladevisual_dict[line[0]] = {"color":pylum_color_dict["Other"], "size":12, "linewidth":0.1, "edgecolor":"#303030"}

tcircle.add_tarc(tarc)
tcircle.set_tarcs()
tcircle.plot_tree(tarc.arc_id, rlim=(0,550), cladevisual_dict=cladevisual_dict, linewidth=0.4, linecolor="#606060")
# plt.show() #or
# tcircle.figure.savefig("tree-1.pdf")

#Set highlights for specific clades
highlight_dict={"Staphylococcus" : {"color":"#808080", "alpha":0.3, "fontsize":5, "label":"Staphylococcus"},
                "Burkholderia"   : {"color":"#B0171F", "alpha":0.3, "fontsize":5, "label":"Burkholderia"},
                "Mycoplasma"     : {"color":"#191970", "alpha":0.3, "fontsize":5, "label":"Mycoplasma"},
                "Synechococcus"  : {"color":"#9ACD32", "alpha":0.3, "fontsize":5, "label":"Synechococcus"}}
tcircle.plot_highlight(tarc.arc_id, highlight_dict)
# plt.show() #or
tcircle.figure.savefig("tree-2.pdf")

#Draw the first rings representing the genetic profiles of ATPases.
from matplotlib.colors import LinearSegmentedColormap
VA_ATP_rings = [{} for _ in range(8)]
F_ATP_rings  = [{} for _ in range(8)]
with open("firstring_attr.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        if int(line[1]) <= 8:
            VA_ATP_rings[int(line[1])-1][line[0]] = 1
        else:
            F_ATP_rings[int(line[1])-9][line[0]]  = 1

cmap = LinearSegmentedColormap.from_list("tmp1", ["#FFFFFF","#dfac03"], N=2)
for i, VA_ATP_ring in enumerate(VA_ATP_rings):
    values = []
    for cladename in tarc.terminal_dict:
        VA_ATP_ring.setdefault(cladename, 0)
        values.append(VA_ATP_ring[cladename])
    tcircle.heatmap(tarc.arc_id, data=values, raxis_range=(550+5*i, 550+5*(i+1)), cmap=cmap)

cmap = LinearSegmentedColormap.from_list("tmp2", ["#FFFFFF","#339933"], N=2)
for i, F_ATP_ring in enumerate(F_ATP_rings):
    values = []
    for cladename in tarc.terminal_dict:
        F_ATP_ring.setdefault(cladename, 0)
        values.append(F_ATP_ring[cladename])
    tcircle.heatmap(tarc.arc_id, data=values, raxis_range=(590+5*i, 590+5*(i+1)), cmap=cmap)

tcircle.figure.savefig("tree-3.pdf")

#Draw the second rings representing the profiles of abundances of fatty acid metabolism.
secondring_dict = {"FA synth init":{}, "FA synth elong":{}, "acyl-CoA synth":{}, "beta-Oxidation":{}, "Ketone biosynth":{}}
with open("secondring_attr.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        secondring_dict[line[2]][line[0]] = float(line[1])

cmap1 = LinearSegmentedColormap.from_list("atp1", ["#FFFFFF","#793a07"], N=500)
cmap2 = LinearSegmentedColormap.from_list("atp2", ["#FFFFFF","#9f1f9f"], N=500)
cmap3 = LinearSegmentedColormap.from_list("atp3", ["#FFFFFF","#0000be"], N=500)
cmap4 = LinearSegmentedColormap.from_list("atp4", ["#FFFFFF","#005500"], N=500)
cmap5 = LinearSegmentedColormap.from_list("atp5", ["#FFFFFF","#b22222"], N=500)
cmaps = [cmap1, cmap2, cmap3, cmap4, cmap5]
for i, key in enumerate(("Ketone biosynth", "beta-Oxidation", "acyl-CoA synth", "FA synth elong", "FA synth init")):
    values = []
    for cladename in tarc.terminal_dict:
        secondring_dict[key].setdefault(cladename, 0)
        values.append(secondring_dict[key][cladename])
    tcircle.heatmap(tarc.arc_id, data=values, raxis_range=(635+40*i, 635+40*(i+1)), cmap=cmaps[i])

tcircle.figure.savefig("tree-4.pdf")

#Draw the genome lengths
clade_genomelength_dict = {}
with open("barplot_attr.csv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(",")
        clade_genomelength_dict[line[0]] = float(line[1])

values     = []
facecolors = []
for cladename in tarc.terminal_dict:
    clade_genomelength_dict.setdefault(cladename, 0)
    values.append(clade_genomelength_dict[cladename])
    facecolors.append(cladevisual_dict[cladename]["color"])
tcircle.barplot(tarc.arc_id, data=values, raxis_range=(840,990), facecolor=facecolors)


tcircle.figure.savefig("tree-5.pdf")


