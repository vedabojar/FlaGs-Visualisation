#!/bin/env python3


# This script generates the FlaGs graphical output as an .html file, the content 
# of which can be downloaded as .svg. The output file also offers some interactivity.


import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import colorsys
import random


#From FlaGs script
def postscriptSize(item):
    if int(item)<1000:
        return(0)
    else:
        return(int(item)/1000)
newQ=0

def operonFamily(item):
    if item==0:
        return ' '
    elif item==center:
        return ' '
    elif item==noProt:
        return ' '
    elif item==noProtP:
        return ' '
    elif item==noColor:
        return ' '
    else:
        return item

# Colours
def random_color(h=None):
    if not h:
        c = int((random.randrange(0,100,5))*3.6)/100
    d = 0.5
    e = 0.5
    return _hls2hex(c, d, e)
        
def _hls2hex(c, d, e):
	return '#%02x%02x%02x' %tuple(map(lambda f: int(f*255),colorsys.hls_to_rgb(c, d, e)))

def outliner (item):
	if item =='#ffffff':
		return '#bebebe'
	elif item =='#f2f2f2':
		return '#008000'
	elif item =='#f2f2f3':
		return '#000080'
	else:
		return item


color={}
colorDict={}

# 1. Merging the two traces into one plot
fig = make_subplots(shared_yaxes = True, shared_xaxes = True)


# 2. Setting the axes. Since the output is a single graph, 
#    the y-axis labels need to be set a specified
data = pd.read_csv('./operon.tsv', sep='\t', header=None, skiprows=None)            
df = pd.DataFrame(data)                                                                                
                                                                                                    
y_axes = df[0].nunique()*40
x_axes = y_axes*5


# 3. All lists
arrowList = []
domainList = []
xList_gene = []
yList_gene = []
xList_domain = []
yList_domain = []
y_tick_marks = []
labels = []
genes = []


# 4. Data file input
main_file = open('operon.tsv','r').read()
eg1 = main_file.split("\n\n\n\n")
y_level_m = 0
for m in eg1:
    if m != '':
        row1 = 0
        entries1 = m.splitlines()
        ndoms=len(entries1)
        y_level_m = y_level_m-10-round(postscriptSize(newQ))
        for entry in entries1:
            items1 = entry.split("\t")
            gene_start = int(items1[5])
            gene_end = int(items1[6])
            gene_length = int(items1[1])
            gene_direction = items1[3]
            dom1_name = int(items1[4])
            accesssion = str(items1[0])
            id1 = str(items1[9][:14])                     

            # 4a. When genes are to small the arrow shape is distorted because the coordinates are too close to each other.
            #     This makes these genes longer to keep the shape of the arrow. 
            if gene_length < 100:
                gene_start = int(items1[5])-50
                gene_end = int(items1[6])+50
            else:
                gene_start = int(items1[5])
                gene_end = int(items1[6])


            #Colours (imported from FlaGs script)
            center=int(dom1_name)+1
            noProt=int(dom1_name)+2
            noProtP=int(dom1_name)+3
            noColor=int(dom1_name)+4
            
            color[noColor]='#ffffff'
            color[center]='#000000'
            color[noProt]='##f2f2f2'
            color[noProtP]='#f2f2f3'


            if dom1_name == 0:
                colorDict[dom1_name]=str('#ffffff')
            elif gene_start == 1:                    
                colorDict[dom1_name]=str('#000000')
            elif 'pseudogene_' in id1:
                colorDict[dom1_name]=str('#f2f2f2') 
            elif 'tRNA_' in id1:
                colorDict[dom1_name]=str('#f2f2f3')
            else:
                if dom1_name not in colorDict:
                    colorDict[dom1_name] = random_color()


            # 4b. Editing the label for each gene in the legend
            protein = (str(id1) + ' ' + '(' + ('Start: {}\tEnd: {}'.format(gene_start, gene_end)) + ')')
            hover_text = 'ID: ' + id1 +'<br>Start: ' + str(gene_start) + '<br>End: ' + str(gene_end)


            # 4c. Drawing the genes as polygons/arrows
            if gene_direction == '-':
                xList_gene = [gene_start+100, gene_start, gene_start+100, gene_end, gene_end, gene_start+100]
                yList_gene = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein)))
                #arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein, legendgroup  = accesssion, legendgrouptitle_text = species)))      # To group the legend
            else:
                xList_gene = [gene_start, gene_start, gene_end-100, gene_end, gene_end-100, gene_start]
                yList_gene = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2] 
                arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein)))
                #arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein, legendgroup  = accesssion, legendgrouptitle_text = species)))      # To group the legend
            

            # 5. Annotating each polygon/arrow with family number
            text_x = gene_start + (gene_length/2.2)
            if dom1_name != 0 and gene_start != 1 and 'pseudogene_' not in id1 and 'tRNA_' not in id1:
                fig.add_annotation(x = text_x, y = y_level_m, xref='x', yref='y', text = dom1_name, font = dict(color = "black", size = 10, family = "Open Sans"), showarrow = False)
            else:
                pass

  
            # 6. Domains file input
            with open ('domains.txt', 'r') as domain_file_org:     # opening domain file
                domain_file = domain_file_org.readlines()[3:-10]     # skipping first 3 rows
                for d in domain_file:
                    if d != '':
                        row2 = 0                                  # row number, goes to y axis
                        y_level_d = y_level_m                     # on what y level to start drawing the arrows
                        entries2 = d.splitlines()                 # 
                        for entry2 in entries2:
                            if entry2 == '':
                                continue # go to end of loop
                            items2 = entry2.split()
                            domain_start = gene_start + ((int(items2[18]))*3)
                            domain_size = ((int(items2[19]))*3)-((int(items2[18]))*3)
                            domain_end = domain_start + domain_size
                            domain_name = str(items2[1])
                            identifier = str(items2[0])
                            id2 = str(items2[4][:14])
                            e_value = items2[7]
                            score = items2[8]
                            accesssion_domains = items2[1]


                            # 6a. Editing the label for each domain in the legend
                            domain = ('     ' + 'Domain: ' + id2 + ' ' + '(' + ('Start: {}\tEnd: {}'.format(domain_start, domain_end)) + ')')


                            # 6b. If a gene has additional information about domains (i.e. same id is found in second file), then these will also be drawn inside the arrow.
                            if id2 == id1:
                                if gene_end-100 < domain_end < gene_end and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, gene_end-100, (gene_end-100) + (gene_end-domain_end), (gene_end-100) + (gene_end-domain_end), gene_end-100, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m+1,  y_level_m-1, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))                                 
                                elif gene_start < domain_start < gene_start+100 and gene_direction == '-':
                                    xList_domain = [gene_start+100, (gene_start+100) - (domain_start-gene_start),  (gene_start+100) - (domain_start-gene_start), gene_start+100, domain_end, domain_end, gene_start+100]
                                    yList_domain = [y_level_m-2, y_level_m-1, y_level_m+1, y_level_m+2,  y_level_m+2, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))                       
                                elif domain_end <= gene_end and domain_start > gene_start+100 and gene_direction == '-':           
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))
                                elif domain_start >= gene_start and gene_direction == '-':
                                    xList_domain = [domain_start+100, domain_start, domain_start+100, domain_end, domain_end, domain_start+100]
                                    yList_domain = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))
                                elif domain_start <= gene_start and domain_end < gene_end-100 and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))
                                elif domain_start <= gene_start and domain_end <= gene_end and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, domain_end-100, domain_end, domain_end-100, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))
                                elif domain_end <= gene_end and domain_start >= gene_start:
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    domainList.append(fig.add_trace(go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = domain)))
    
                                #print(accesssion_operons + "(" + id1 + ")" + " : " + id2)
                                print(identifier)


            # 7. Setting the y labels i.e. the organism name and accession nr etc.
            y_tick_marks += [y_level_m]
            labels += [items1[0]]
            genes += [items1[9]]

            row1 = row1+1


# 8. Changing the download format of the .html as .svg instead of the defaul .png
config = {'toImageButtonOptions': {'format': 'svg','filename': 'FlaGs','scale': 1}}


# 9. Graph layout
fig.update_xaxes(visible = False)
fig.update_yaxes(visible = True, showgrid = False, showline = False, autorange = True, automargin = True, showticklabels = True, tickvals = y_tick_marks, ticktext = labels, ticklen = 20, tickmode = 'array', titlefont = dict(family = 'Open Sans', size = 8))
fig.update_layout(autosize=False, width=x_axes, height=y_axes, margin=dict(l=100,r=500,b=100,t=100,pad=100),paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', showlegend = True)

fig.show(config=config)

#fig.write_html("test.html")
