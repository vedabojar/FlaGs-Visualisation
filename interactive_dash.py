#!/bin/env python3


# This script generates the FlaGs graphical output with imporved interactivity 
# where the user can filter, select, click, and get more info about the data. 
#
# Run the applicatuib with "python app.py" and visit 
# http://127.0.0.1:8050/ in your web browser.


import dash_bio as dashbio
from matplotlib.pyplot import plot, show, yticks
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from dash import html, dcc, dash_table
import pandas as pd
import colorsys
import random
import re
import dash_daq as daq
import dash_bootstrap_components as dbc
import dash
from dash.dependencies import Input, Output, State
from collections import Counter  
import plotly.figure_factory as ff
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline



#app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)

server = app.server

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
fig1 = make_subplots(shared_yaxes = True, shared_xaxes = True)
fig2 = make_subplots(shared_yaxes = True, shared_xaxes = True)


# 3. All lists
arrowList = []
domainList = []
xList_gene = []
yList_gene = []
xList_domain = []
yList_domain = []
xList_domain_line = []
yList_domain_line = []
y_tick_marks = []
labels = []
gene_list = []
domain_list = []
id1 = []
specie_name = []
annotations_list = []
trace1A_list = []   # traces for the operons - operon graph
trace1B_list = []   # traces for the domains - operon grpah
traceB_list = []    # traces for the domains - domain graph
traceB_map = dict()
accession_operons_List = []
accession_domains_List = []

# 4. Data file input
main_file = open('dm_operon_short.tsv','r').read()
eg1 = main_file.split("\n\n\n\n")
y_level_m = 0
for m in eg1:
    if m != '':
        row1 = 0
        entries1 = m.splitlines()
        ndoms = len(entries1)
        y_level_m = y_level_m-10-round(postscriptSize(newQ))
        for entry in entries1:
            entry = re.sub("\s\s+"," ", entry)
            entry = entry.replace(" ", "\t")
            items1 = entry.split("\t")
            gene_start = int(items1[5])
            gene_end = int(items1[6]) 
            gene_length = int(items1[1])
            gene_direction = items1[3]
            dom1_name = int(items1[4])
            id1 = str(items1[9][:14]) 
            accesssion_operons = str(items1[0])
            accession_operons_List.append(accesssion_operons)                 

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
            hover_text = 'ID: ' + id1 


            # 4c. Drawing the genes as polygons/arrows
            if gene_direction == '-':
                xList_gene = [gene_start+100, gene_start, gene_start+100, gene_end, gene_end, gene_start+100]
                yList_gene = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                #arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein), row = 1, col = 1))
                trace1A_list.append((go.Scatter(x = xList_gene, y = yList_gene, fill="toself", hovertext = hover_text, fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), hoveron = "fills", name = id1, mode = 'lines+text'), len(fig1.data)))
                fig1.add_trace(trace1A_list[-1][0], row = 1, col = 1)
            else:
                xList_gene = [gene_start, gene_start, gene_end-100, gene_end, gene_end-100, gene_start]
                yList_gene = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2] 
                #arrowList.append(fig.add_trace(go.Scatter(x = xList_gene, y = yList_gene, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), mode = 'lines+text', name = protein), row = 1, col = 1))
                trace1A_list.append((go.Scatter(x = xList_gene, y = yList_gene, fill="toself", hovertext = hover_text, fillcolor=colorDict[dom1_name], opacity = 0.5, line=dict(color = outliner(colorDict[dom1_name])), hoveron = "fills", name = id1, mode = 'lines+text'), len(fig1.data)))
                fig1.add_trace(trace1A_list[-1][0], row = 1, col = 1)
            

            # 5. Annotating each polygon/arrow with family number
            text_x = gene_start + (gene_length/2.2)
            if dom1_name != 0 and gene_start != 1 and 'pseudogene_' not in id1 and 'RNA_' not in id1:
                annotations_list.append(fig1.add_annotation(x = text_x, y = y_level_m, xref='x', yref='y', text = dom1_name, font = dict(color = "black", size = 8, family = "Open Sans"), showarrow = False))
            else:
                pass
            
  
            # 6. Domains file input
            domainList = []
            coordList = []
            scoreList = []

        
            data_2 = pd.read_csv('./dm_output.txt', sep = "\s+|\t+|\s+\t+|\t+\s+", header = 1,  usecols=range(24), engine='python')      #NOTE! Table might be incomplete,       
            df_2 = pd.DataFrame(data_2)       

            with open ('dm_output.txt', 'r') as domain_file_org:     # opening domain file
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
                            accession_domains_List.append(accesssion_domains)
                            domainList.append(domain_name)
                            
                            # 6a. Editing the label for each domain in the legend
                            domain = ('     ' + 'Domain: ' + id2 + ' ' + '(' + ('Start: {}\tEnd: {}'.format(domain_start, domain_end)) + ')')

                            # 6b. If a gene has additional information about domains (i.e. same id is found in second file), then these will also be drawn inside the arrow.
                            if id2 == id1:
                                if gene_end-100 < domain_end < gene_end and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, gene_end-100, (gene_end-100) + (gene_end-domain_end), (gene_end-100) + (gene_end-domain_end), gene_end-100, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m+1,  y_level_m-1, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)                                       
                                elif gene_start < domain_start < gene_start+100 and gene_direction == '-':
                                    xList_domain = [gene_start+100, (gene_start+100) - (domain_start-gene_start),  (gene_start+100) - (domain_start-gene_start), gene_start+100, domain_end, domain_end, gene_start+100]
                                    yList_domain = [y_level_m-2, y_level_m-1, y_level_m+1, y_level_m+2,  y_level_m+2, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))                                       
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)                            
                                elif domain_end <= gene_end and domain_start > gene_start+100 and gene_direction == '-':           
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)
                                elif domain_start >= gene_start and gene_direction == '-':
                                    xList_domain = [domain_start+100, domain_start, domain_start+100, domain_end, domain_end, domain_start+100]
                                    yList_domain = [y_level_m-2, y_level_m, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))
                                elif domain_start <= gene_start and domain_end < gene_end-100 and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))                                          
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)
                                elif domain_start <= gene_start and domain_end <= gene_end and gene_direction == '+':
                                    xList_domain = [domain_start, domain_start, domain_end-100, domain_end, domain_end-100, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))                                       
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)
                                elif domain_end <= gene_end and domain_start >= gene_start:
                                    xList_domain = [domain_start, domain_start, domain_end, domain_end, domain_start]
                                    yList_domain = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                                    trace1B_list.append((go.Scatter(x=xList_domain, y=yList_domain, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', text = domain_name, name = identifier), len(fig1.data)))        
                                    fig1.add_trace(trace1B_list[-1][0], row = 1, col = 1)
    
                                #print(accesssion_operons + "(" + id1 + ")" + " : " + id2)
                                print(identifier)


                        # # 6c. If a gene has additional information about domains (i.e. same id is found in second file), then these will also be drawn inside the arrow.
                        # if id1 == id2:
                        #     if domain_end != gene_end and domain_start != gene_start:
                        #         xList_domain_line = [domain_start, domain_start, domain_end, domain_end, domain_start]
                        #         yList_domain_line = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                        #         traceB_list.append((go.Scatter(x=xList_domain_line, y=yList_domain_line, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = id2), len(fig2.data)))     
                        #         fig2.add_trace(traceB_list[-1][0], row = 1, col = 1)
                        #     elif domain_end == gene_end and gene_direction == '-':
                        #         xList_domain_line = [domain_start, domain_start, domain_end, domain_end, domain_start]
                        #         yList_domain_line = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                        #         traceB_list.append((go.Scatter(x=xList_domain_line, y=yList_domain_line, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = id2), len(fig2.data)))
                        #         fig2.add_trace(traceB_list[-1][0], row = 1, col = 1)
                        #     elif domain_start == gene_start and gene_direction == '+':
                        #         xList_domain_line = [domain_start, domain_start, domain_end, domain_end, domain_start]
                        #         yList_domain_line = [y_level_m-2, y_level_m+2, y_level_m+2, y_level_m-2, y_level_m-2]
                        #         traceB_list.append((go.Scatter(x=xList_domain_line, y=yList_domain_line, fill="toself", hoverinfo = 'none', fillcolor=colorDict[dom1_name], line=dict(color=colorDict[dom1_name]), opacity = 0.5, mode='lines', name = id2), len(fig2.data)))
                        #         fig2.add_trace(traceB_list[-1][0], row = 1, col = 1)                                                  

                                # Group all fig.2 traces that have the same ID in column = 0  
                                if id2 not in traceB_map.keys():
                                    traceB_map[id2] = []
                                traceB_map[id2].append(len(fig2.data) - 1)

                                # 6c. Annotating each domain with name
                                text_x = domain_start + (domain_size/2)
                                fig2.add_annotation(x = text_x, y = y_level_d+5, xref='x', yref='y', text = domain_name, font = dict(color = "black", size = 8, family = "Open Sans"), showarrow = False)


            # 7. Setting the y labels i.e. the organism name and accession nr etc.
            y_tick_marks += [y_level_m]
            labels += [items1[0]]
            gene_list += [items1[9][:14]]


            row1 = row1+1


#Alignment
from Bio.Align.Applications import MafftCommandline
in_file = "sequences.fasta"
mafft_cline = MafftCommandline(input=in_file)
print(mafft_cline)
stdout, stderr = mafft_cline()
with open("aligned.fasta", "w") as handle:
    handle.write(stdout)
from Bio import AlignIO
aligned = AlignIO.read("aligned.fasta", "fasta")
alignment_data = open('aligned.fasta', 'r').read()
print(alignment_data)
#output.write(stdout) 

# ------------------- DOMAIN GRAPH -------------------

# 1. Protein file

protein_file = pd.read_csv("./dm_operon_short.tsv", sep = "\t", header = None, engine='python')

protein_dict = {
    "accession": [], 
    "Task": [], 
    "Start": [],
    "Finish": [],
    "E-value": [],
    "Score": [],
}

for str in protein_file.iloc[:, 0]:
    protein_dict["Task"].append(str.split('#')[0])

for str in protein_file.iloc[:, 9]:
    protein_dict["accession"].append(str.split('#')[0])
    
    
protein_dict["Start"] = [0]*len(protein_file)
protein_dict["Finish"] = list(protein_file.iloc[:, 6].astype(int)/3 - protein_file.iloc[:, 5].astype(int)/3)
protein_dict["E-value"] = [-1000000]*len(protein_file)
protein_dict["Score"] = [1000000000]*len(protein_file)

protein_df = pd.DataFrame(protein_dict)



# 2. Domain file 

count = 0
newFileStr = ""
with open("dm_output.txt") as fp:
    Lines = fp.readlines()
    for line in Lines:
        count += 1
        newFileStr += re.sub('[\t, +]+', ',', line)           # Crazy separators need relpacing
        
text_file = open("dm_output.csv", "w")
n = text_file.write(newFileStr)
text_file.close()

domain_file = pd.read_csv('./dm_output.csv', usecols=range(24), header = 1, comment='#', engine='python')        # Even with replacement of sep, the last column
                                                                                                # is read wrong, therefore columns up to 23 are read
domain_dict = { 
    "accession": [], 
    "Task": [], 
    "Start": [],
    "Finish": [],
    "E-value": [],
    "Score": [],
    "identifier": []
}

for str in domain_file.iloc[:, 4]:
    domain_dict["accession"].append(str.split('|')[1])   

domain_dict["Start"] = list(domain_file.iloc[:, 20])
domain_dict["Finish"] = list(domain_file.iloc[:, 21])
domain_dict["Task"] = list(domain_file.iloc[:, 1].replace('_', ' ', regex=True))            
domain_dict["E-value"] = list(domain_file.iloc[:, 7])
domain_dict["Score"] = list(domain_file.iloc[:, 8])
domain_dict["identifier"] = list(domain_file.iloc[:, 0])

domain_df = pd.DataFrame(domain_dict)





datatable_color = ['#B0171F',       # NOTE! MUST ADD COLORS 
        '#DC143C',
        '#FFB6C1',
        '#FFAEB9',
        '#EEA2AD',
        '#CD8C95',
        '#8B5F65',
        '#FFC0CB',
        '#FFB5C5',
        '#EEA9B8',
        '#CD919E',
        '#8B636C',
        '#DB7093',
        '#FF82AB',
        '#EE799F',
        '#CD6889',
        '#8B475D',
        '#FFF0F5',
        '#EEE0E5',
        '#CDC1C5',
        '#8B8386',
        '#FF3E96',
        '#EE3A8C',
        '#CD3278',
        '#8B2252',
        '#FF69B4',
        '#FF6EB4',
        '#EE6AA7',
        '#CD6090',
        '#8B3A62',
        '#872657',
        '#FF1493',
        '#EE1289',
        '#CD1076',
        '#8B0A50',
        '#FF34B3',
        '#EE30A7',
        '#CD2990',
        '#8B1C62',
        '#C71585',
        '#D02090',
        '#DA70D6',
        '#FF83FA']


# 8. Changing the download format of the .html as .svg instead of the defaul .png
config = {'toImageButtonOptions': {'format': 'svg','filename': 'FlaGs','scale': 1}}


y_axes_operon = ((len(Counter(accession_operons_List).keys()))*2.0)      
x_axes_operon = y_axes_operon * 2.0     

y_axes_domain = ((len(Counter(accession_domains_List).keys()))*2.0)      
x_axes_domain = y_axes_domain * 2.0 


# 9A. Operon graph layout
fig1.update_xaxes(visible = False)
fig1.update_yaxes(visible = True, 
                showgrid = False, 
                showline = False, 
                autorange = True, 
                automargin = True, 
                showticklabels = True, 
                tickvals = y_tick_marks, 
                ticktext = labels, 
                ticklen = 20, 
                tickmode = 'array', 
                tickfont = dict(
                            family = 'arial', 
                            size = 10))
fig1.update_layout(autosize=False, 
                width = x_axes_operon, 
                height = y_axes_operon, 
                margin = dict(l=10,r=10,b=10,t=5,pad=10),
                paper_bgcolor = '#ffffff', 
                plot_bgcolor = 'rgba(0,0,0,0)', 
                showlegend = False, 
                clickmode = 'event')
fig1.update_traces(marker_size = 20)


# 9B. Domain graph layout
fig2.update_xaxes(visible = False)
fig2.update_yaxes(visible = False, 
                showgrid = False, 
                showline = False, 
                autorange = True, 
                showticklabels = False)
fig2.update_layout(autosize=True, 
                #width=x_axes, 
                #height=100, 
                margin=dict(l=10,r=10,b=10,t=5,pad=10),
                paper_bgcolor='#ffffff', 
                plot_bgcolor='rgba(0,0,0,0)', 
                showlegend = False, 
                clickmode = 'event+select')
fig2.update_traces(marker_size = 20)


#fig.show(config=config)
#fig.write_html("test.html")

DEFAULT_COLORSCALE = ["#2a4858", "#265465", "#1e6172", "#106e7c", "#007b84", \
	"#00898a", "#00968e", "#19a390", "#31b08f", "#4abd8c", "#64c988", \
	"#80d482", "#9cdf7c", "#bae976", "#d9f271", "#fafa6e"]

DEFAULT_OPACITY = 0.8

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}



app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            dbc.Row([
                html.H3('FlaGs Interactive Output', 
                        style={'margin-left': '40px',
                            'margin-bottom': '60px',
                            'margin-top' : '50px'}),
            ]),
            dbc.Row([
                html.Details([
                    html.Summary('Domain Plot'),
                        html.Label(['Sort by:'],
                                style={'font-weight': 'bold', 
                                    'margin-left' : '20px',
                                    'margin-top' : '20px', 
                                    'backgroundColor':'#f9f9f9'}),
                        dcc.RadioItems(
                                id='radio',
                                options=[
                                        {'label': ' Sequence order', 'value': 'Sequence order'},
                                        {'label': ' Lowest E-value', 'value': 'Lowest E-value'},
                                ],
                                value='Sequence order',
                                labelStyle = {'display': 'inline-block',
                                            'margin-right': 20},
                                style={"width": "60%",
                                    'margin-left' : '20px',
                                    'margin-top' : '10px',
                                    'margin-right' : '10px',
                                    'backgroundColor':'#f9f9f9'}),
                        html.Hr(style={"width": "95%",
                                    'margin-left' : '20px'}),                   
                        dcc.Graph(
                                id = 'domain-plot',
                                figure = {},
                                style = {'display': 'block',
                                        'margin-left' : '20px',
                                        'margin-right' : '5px',
                                        'margin-top' : '20px',
                                        'margin-bottom' : '10px',
                                        'width' : '95%',
                                        'height' : "auto"}),
                        ])
            ], style={'backgroundColor':'#f9f9f9', 
                    'margin-left' : '40px', 
                    'margin-top' : '10px',
                    'position': 'sticky',
                    'top': 0,
                    'z-index': 999}),
            dbc.Row([
                html.H2(id = 'clickdata-output'),
                # html.Tr([
                #     html.Td("Hide/Show Domains", 
                #             style={'display': 'inline-block',
                #                     'margin-left' : '20px',
                #                     'margin-right' : '10px',
                #                     'width' : '95%',
                #                     'height' : '40px',
                #                     'vertical-align' : 'middle'}),       #NOTE: Align text in center
                #     html.Td(
                #         daq.BooleanSwitch(
                #             id = 'Hide/Show Domains switch', 
                #             on = True, 
                #             color = "#2a4858",
                #             style = {'display': 'inline-block',
                #                     'margin-left' : '30px',
                #                     'margin-right' : '10px',
                #                     'margin-top' : '10px',
                #                     'width' : '95%',
                #                     'height' : '40px',
                #                     'postion': 'relative', 
                #                     'top': 0,
                #                     'z-index': 0}
                #                 )
                #             )
                #         ]),
                # html.Hr(style={"width": "93%",
                #             'margin-left' : '30px',
                #             'postion': 'relative', 
                #             'top': 0,
                #             'z-index': 0}), 
                dcc.Graph(
                    id = 'operon-plot',
                    animate = False,
                    figure = fig1,
                    responsive = True,
                    style = {'display': 'block',
                            'margin-left' : '20px',
                            'margin-right' : '5px',
                            'margin-top' : '5px',
                            'margin-bottom' : '10px',
                            'width' : '95%',
                            'height' : '300px'}),
            ], style={'backgroundColor':'#f9f9f9', 
                    'margin-left' : '40px', 
                    'margin-top' : '10px',
                    'margin-bottom' : '10px',
                    'postion': 'relative', 
                    'top': 0,
                    'z-index': 0}),
            dbc.Row([
                html.Details([
                    html.Summary('Alignment Chart'),       
                        html.Div([
                                dashbio.AlignmentChart(
                                    id='my-default-alignment-viewer',
                                    data=alignment_data,
                                    height=500,
                                    tilewidth=40,
                                ),
                        ])
                ]),
            ], style={'backgroundColor':'#f9f9f9', 
                    'margin-left' : '40px', 
                    'margin-top' : '10px', 
                    'margin-bottom' : '10px'}),
        ], width = 8),
        dbc.Col([
            dbc.Row([
                html.Div([
                    dcc.Dropdown(
                        id = 'drop-down-proteins',
                        placeholder = 'Select a protein for more details',
                        value = gene_list,
                        options = gene_list,
                        #color = '#99B2B8',
                        style = {'display': 'inline-block',
                                'margin-left' : '2px',
                                'margin-right' : '3px',
                                'margin-top' : '5px',
                                'width' : '100%',
                                'height' : '40px'}),
                    html.Br(),  
                    #dbc.Alert(
                        # id = "Gene info",
                        # children = "Select/Click on a gene to see its details here",
                        # color = "secondary", 
                        # style = {'display': 'inline-block',
                        #         'margin-left' : '5px',
                        #         'margin-top' : '10px',
                        #         'width' : '95%'}),
                    ])
                ], "mb-5"),
            dbc.Row([
                dash_table.DataTable(
                    id = 'domain-checkbox-table',
                    row_selectable = 'multi',
                    selected_rows = [],
                    editable = True,
                    style_cell_conditional=[{'if': {'column_id': c},'textAlign': 'left'} for c in ["Task"]],
                    #filter_action = 'native',
                    sort_action = 'native',
                    style_cell = {'fontSize': 9, 
                                'font-family':'arial',
                                'minWidth': '60px',   # 92px
                                'width': '60px', 
                                'maxWidth': '60px',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis'}, 
                    style_data = {'whiteSpace': 'normal', 
                                'height': 'auto', 
                                'width': '100%'},
                    style_table = {'margin-left' : '2px',
                                'width' : '100%'},            
                    fill_width = False), 
                ], className="h-auto"),
            html.Hr(style={"width": "93%",
                        'margin-left' : '10px'}), 
            html.Tr([
                    html.Td("Hide/Show Selected Domain", 
                            style={'display': 'inline-block',
                                    'margin-left' : '10px',
                                    'margin-right' : '10px',
                                    'width' : '95%',
                                    'height' : '40px',
                                    'vertical-align' : 'middle', 
                                    'fontSize': 12}),       #NOTE: Align text in center
                    html.Td(
                        daq.BooleanSwitch(
                            id = "Hide/Show Selected Domain", 
                            on = True, 
                            color = "#2a4858",
                            style = {'display': 'inline-block',
                                    'margin-left' : '20px',
                                    'margin-right' : '10px',
                                    'margin-top' : '10px',
                                    'width' : '95%',
                                    'height' : '40px'}))
                        ]),
            html.Tr([
                    html.Td("Hide/Show Selected Domain Annotations", 
                            style={'display': 'inline-block',
                                    'margin-left' : '10px',
                                    'margin-right' : '10px',
                                    'width' : '95%',
                                    'height' : '40px',
                                    'vertical-align' : 'middle',
                                    'fontSize': 12}),       #NOTE: Align text in center
                    html.Td(
                        daq.BooleanSwitch(
                            id = 'Hide/Show Selected Domain Annotations', 
                            on = True, 
                            color = "#2a4858",
                            style = {'display': 'inline-block',
                                    'margin-left' : '20px',
                                    'margin-right' : '10px',
                                    'margin-top' : '10px',
                                    'width' : '95%',
                                    'height' : '40px'}))
                        ]),
        ], style={'backgroundColor':'#f9f9f9', 
                    'margin-left' : '40px',
                    'margin-top' : '150px',
                    'margin-bottom' : '10px'},
        width = 3),
    ])
], fluid = True, style={'backgroundColor':'#f2f2f1'})








# HEATMAP
@app.callback(
    Output('default-alignment-viewer-output', 'children'),
    Input('my-default-alignment-viewer', 'eventDatum')
)
def update_output(value):
    if value is None:
        return 'No data.'
    return str(value)


# 2. Toggle Hide/Show switches:
#2A. Domains
@app.callback(
    [   Output('operon-plot', 'figure'),
        Output('drop-down-proteins', 'value')
    ],
    #Input('Hide/Show Domains switch', 'on'),)
    #Input('domain-checkbox-table', 'derived_virtual_data'),
    [   Input('operon-plot', 'clickData'),
    ]
)
def update_graph(clicked_protein):    #index 
    # if show_hide == True:
    #     for t in range(len(trace1B_list)):
    #         fig1.data[trace1B_list[t][1]]["opacity"] = 0.5
    #     return fig1.update_traces()
    # else:
    #     for t in range(len(trace1B_list)):
    #         fig1.data[trace1B_list[t][1]]["opacity"] = 0.0
    #     return fig1.update_traces()
    if clicked_protein is not None:
        #for t in range(len(trace1A_list)):
            curve_number = clicked_protein["points"][0]["curveNumber"]
            trace_name = fig1['data'][curve_number]['name']
            fig1['data'][curve_number]['line']['color'] = '#000000'
    return fig1.update_traces(), trace_name

    if index is not None:
        for i in index:
            id_table = i["identifier"]
            #print(id_table)
            t = 0
            for t in range(len(trace1B_list)): 
                id_graph = int(trace1B_list[t][0]["name"])
                #domain_tag = int(trace1B_list[t][0]["text"])
                #trace_opacity = int(trace1B_list[t][0]["opacity"])
                if id_table == id_graph:
                    return print(id_graph)
                    #print(trace_opacity)
                    #fig1.data[trace1B_list[t][0]]["opacity"] = 0.5
                    #return fig1.update_traces()
                else: 
                    t = t + 1

            # for t in range(len(trace1B_list)): 
            #     for name_2 in range(len(trace1B_list[t][0]["text"])):
            #         if name == name_2:
            #             fig1.data[trace1B_list[t][0]]["opacity"] = 0.0
            #return fig1.update_trace()



    # if prot is not None:
    #     for prot in range(len(domain_df)):


    # if checked_domain is not None:
    #     for t in domain_df["Task"]:
    #         if t == checked_domain: 
    #             fig1.data[trace1B_list[t][1]]["opacity"] = 0.0
    #         else:
    #             fig1.data[trace1B_list[t][1]]["opacity"] = 0.5
    #     return fig1.update_traces()

    # for mark_protein in range(len(trace1A_list)):
    #     fig1.data[trace1A_list[mark_protein][0]["line"]["color"]] = "black"
    # return fig1.update_traces()

    

        

# 2A. 2nd Graph
@app.callback(
    Output('domain-plot', 'figure'),
    [   Input('drop-down-proteins', 'value'), 
        Input('radio', 'value'),
    ]
)
def domain_fig(selected_protein, sort):
    for prot in protein_df["accession"]:
        if prot == selected_protein:
            merged_df = pd.concat([protein_df.loc[protein_df['accession'] == prot], domain_df.loc[domain_df['accession'] == prot]], axis = 0)
            if sort == 'Lowest E-value':
                plot_df = merged_df.sort_values(by=['E-value'], ascending = True)
            else:
                sort == 'Sequence order'
                plot_df = merged_df.sort_values(by=['Start'], ascending = True)

            y = len(plot_df['Task'].unique())
            title = ('Protein: {}'.format(prot))
            fig_g = ff.create_gantt(plot_df, colors = datatable_color, index_col='Task', bar_width = 0.2, height = y * 40)
            fig_g.layout.title = title
            fig_g.update_layout(paper_bgcolor='#ffffff', 
                                plot_bgcolor='rgba(0,0,0,0)', 
                                margin=dict(l=10,r=10,b=5,t=35,pad=10),
                                xaxis_type='linear',
                                xaxis = dict(visible = False),
                                yaxis = dict(visible = True, 
                                            autorange = 'reversed', 
                                            showticklabels = True, 
                                            automargin = True))
    return fig_g                



#Display in table
@app.callback(
    Output('domain-checkbox-table', 'data'),
    Output('domain-checkbox-table', 'columns'),
    Input('drop-down-proteins', 'value')) 
    #Input('operon-plot', 'clickData'))
def update_table(selected_protein):
    #clicked_protein == 
    for prot in protein_df["accession"]:
        if prot == selected_protein:
            table_df = domain_df[domain_df['accession'] == prot]
            columns=[{'name': i, 'id': i, "selectable": True} for i in table_df.loc[:,("Task", "E-value", "Score", "Start", "Finish")]]
            data = table_df.to_dict('records')
            return data, columns
            

# @app.callback(
#     Output('operon-plot','figure'),
#     Input('domain-checkbox-table', 'selected_rows'))
# def update_domains(selected_protein):
#     for prot in protein_df["accession"]:
#         for t in range(len(trace1B_list)):
#             if prot == 



# ---------- Clickable Traces --------------
# @app.callback(
#     Output('drop-down-proteins', 'value'),
#     Input('operon-plot', 'clickData'))
# def click_on_protein(clicked_protein):
#     if clicked_protein is not None:
#         curve_number = clicked_protein["points"][0]["curveNumber"]
#         trace_name = fig1['data'][curve_number]['name']
#     return trace_name



# Select gene from list 
@app.callback(
    Output('Gene info', 'children'),
    Input('drop-down-proteins', 'value'))
def display_info(selected_gene):
    if type(selected_gene) != str:
        return
    for gene, id in trace1A_list:
        if gene["name"] == selected_gene:
            gene_info_text = []
            gene_info_text.append(html.P("ID: " + fig1.data[id]["name"][:16]))
            gene_info_text.append(html.P("Start: " + str(fig1.data[id]['x'][0])))
            gene_info_text.append(html.P("Stop: " + str(fig1.data[id]['x'][2])))
    return gene_info_text 





if __name__ == '__main__':
    app.run_server(host="127.0.0.1", port="8050", debug=True)
