#!/bin/env python3


import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, tight_layout          # NOTE! Add this
import pandas as pd
import colorsys
import random
from collections import Counter                             # NOTE! Add this
1
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
arrowList = []
#domainList = []        # Needed if step 4 (showing domains) is used
accession_List = []

# 2. Drawing the plots
fig, ax = plt.subplots(1, 2, sharey = 'row', gridspec_kw={'width_ratios': [5, 15]})


# 3. Operon (main) file input 
main_file = open('operon.tsv','r').read()
eg1 = main_file.split("\n\n\n\n")
y_level_m = 0
for m in eg1:
    if m != '':
        row = 0
        entries1 = m.splitlines()
        ndoms=len(entries1)
        y_level_m = y_level_m-7-round(postscriptSize(newQ))
        for entry in entries1:
            items1 = entry.split("\t")
            x_gene_start = int(items1[5])
            x_gene_end = int(items1[6])
            dx_gene_length = int(items1[1])
            gene_direction = items1[3]
            dom1_name = int(items1[4])
            id1 = str(items1[9])
            accession = str(items1[0])                      # NOTE! Add this
            accession_List.append(accession)                # NOTE! Add this         

            # 3a. When genes are to small the arrow shape is distorted because the coordinates are too close to each other.
            #     This makes these genes longer to keep the shape of the arrow. 
            if dx_gene_length < 80:
                dx_gene_length = int(items1[1])*2.5
            else:
                dx_gene_length = int(items1[1])

            # 3b. Gene diretion and the pointing the arrow head in the correct direction
            if gene_direction == '-':
                x_gene_start = int(items1[6])
                x_gene_end = int(items1[5])
                dx_gene_length = dx_gene_length*(-1)
            else:
                x_gene_start = int(items1[5])
                x_gene_end = int(items1[6])
                dx_gene_length = dx_gene_length


            # Colours (imported from FlaGs script)
            center=int(dom1_name)+1
            noProt=int(dom1_name)+2
            noProtP=int(dom1_name)+3
            noColor=int(dom1_name)+4
            
            color[noColor]='#ffffff'
            color[center]='#000000'
            color[noProt]='#f2f2f2'
            color[noProtP]='#f2f2f3'


            if dom1_name == 0:
                colorDict[dom1_name]=str('#ffffff')
            elif x_gene_start == 1:                    
                colorDict[dom1_name]=str('#000000')
            elif 'pseudogene_' in id1:
                colorDict[dom1_name]=str('#f2f2f2') 
            elif 'tRNA_' in id1:
                colorDict[dom1_name]=str('#f2f2f3')
            else:
                if dom1_name not in colorDict:
                    colorDict[dom1_name] = random_color()



            # 3c. Drawing the genes as arrows. 
            arrowList.append(ax[1].arrow(x=x_gene_start, y=y_level_m, dx=dx_gene_length, dy=0, width=5, head_width=5, length_includes_head = True, head_length = 120, facecolor = colorDict[dom1_name], edgecolor = outliner(colorDict[dom1_name]), alpha=1))

            # # 4 Domains file input (This step will be skiped?)
            # domain_file = open('test2_doms.txt', 'r').read()
            # eg2 = domain_file.split("\n\n\n\n")
            # id1 = str(items1[9])
            # #while != '':
            # for d in eg2:
            #     if d != '':
            #         row1 = 0
            #         entries2=d.splitlines()
            #         ndoms=len(entries2)
            #         y_level_d = y_level_m 
            #         for entry2 in entries2:
            #             if entry2 == '':
            #                 continue # go to end of loop
            #             items2 = entry2.split('\t')
            #             x_domain_start = int(items2[3])
            #             x_domain_end = int(items2[4])
            #             dx_domain_size = int(x_domain_end)-int(x_domain_start)
            #             id2 = str(items2[0])
            #             domain_name = str(items2[2])
    
            #             # 4a. If there 
            #             if id1 == id2:
            #                 if x_domain_end != x_gene_end:
            #                     domainList.append(ax[1].arrow(x=x_domain_start, y=y_level_d, dx=dx_domain_size, dy=0, width=5, head_width = 5, length_includes_head=True, head_length=0, facecolor=colorDict[dom1_name], edgecolor = colorDict[dom1_name], alpha=0.4))
            #                 else:
            #                     domainList.append(ax[1].arrow(x=x_domain_start, y=y_level_d, dx=dx_domain_size, dy=0, width=5, head_width = 5, length_includes_head=True, head_length=120, facecolor=colorDict[dom1_name], edgecolor = colorDict[dom1_name], alpha=0.4))
            #             else:
            #                 pass   

            # 5. Adding the family number inside the gene/arrow
            text_x = x_gene_start + (dx_gene_length/2)
 
            if dom1_name != 0 and x_gene_start != 1 and 'pseudogene_' not in id1 and 'tRNA_' not in id1:
                ax[1].text(text_x, y_level_m-1, s = dom1_name, horizontalalignment='center', color = 'k', fontweight = 'bold', font = {'family' : 'Helvetica','size'   : 8})
            else:
                pass

            # 6. Adding the second plot (left-hand side) with the organism name and accesion nr etc. 
            ptnstats = entries1[0].split("\t")
            org = ptnstats[0][:ptnstats[0].index('|')]+ptnstats[0][ptnstats[0].index('|'):].replace('_',' ')
            ax[0].text(0.7, y_level_m, org, horizontalalignment='center', color = '#000000', font = {'family' : 'Helvetica','size':9})
            ax[0].set_axis_off()
            
            row = row+1


#NOTE! You can copy paste these lines
numb_queries = len(Counter(accession_List).keys())

if numb_queries <= 10:
    y_new = (numb_queries * 0.5)
    x_new = y_new + 10.0
elif 11 < numb_queries <= 499:
    y_new = (numb_queries * 0.02)
    x_new = y_new + 10.0
elif numb_queries > 500:
    y_new = (numb_queries * 0.01)
    x_new = y_new + 10.0

fig.set_figheight(y_new)
fig.set_figwidth(x_new)
plt.tight_layout()
plt.axis('off')
plt.show()

#fig.savefig('test.svg', bbox_inches='tight')




# ATTEMPT TO MAKE THE OUTPUT INTERACTIVE


# 1. Hover to arrow highlight

#def hover(event):
#    if event.inaxes == ax[1]:
#        for i, a in enumerate(domainList):
#            if a.contains_point([event.x, event.y]):
#                a.set_alpha(1.0)
#            else:
#                a.set_alpha(0.5)
#        fig.canvas.draw_idle()


#fig.canvas.mpl_connect("motion_notify_event", hover)

#annot = ax[1].annotate("", xy=(0,0), xytext=([40,40]), textcoords='offset points', bbox=dict(boxstyle='square', fc='gainsboro', ec='gainsboro', alpha=0.7), arrowprops=dict(arrowstyle='wedge', connectionstyle='arc3', fc='gainsboro', edgecolor = 'gainsboro', alpha = 0.7))
#annot.set_visible(False)




# 2. Click for more information

#coord = []
#def on_click (event):
#    global coord 
#    coord.append((event.xdata, event.ydata))
#    x = event.xdata
#    y = event.ydata
#    annot.xy = (x,y)
#    text_click = ('ID: {}\nDomain: {}\nStart: {}\nStop: {}'.format((items2[0]), (items2[2]), items2[3], items2[4]))
#    annot.set_text(text_click)
#    annot.set_visible(True)
#    fig.canvas.draw()
#fig.canvas.mpl_connect("button_press_event", on_click)

#plt.subplots_adjust(left=0.25)

# Checkbutton widget
#labels = (domainList, arrowList)
#activated = [True, False, False]
#axCheckButton = plt.axes([0.03, 0.4, 0.15, 0.15])
#checkbox = CheckButtons(axCheckButton, labels, activated)


#ax[1].annotate("Domain A",size=20, va="center", ha="center", bbox=dict(boxstyle="square", fc="w", alpha = 0.5),arrowprops=dict(arrowstyle="wedge",connectionstyle="arc3",fc="w", alpha=0.5),)


#fig.canvas.mpl_connect('add', show_annotation)

#cr = mplcursors.cursor(hover=True, highlight=True)
#@cr.connect("add")
#def on_add (sel):
#    sel.annotation.set_text('Domain\nID: {}\nLength: {}'.format(df.loc[sel.target.index,9], df.loc[sel.target.index,1]))
#    sel.annotation.get_bbox_patch().set(boxstyle="square", fc="white", alpha=0.5)
    #sel.annotation.get_arrowstyle().set(arrowstyle="wedge", conncectionstyle="angle3")
    #sel.highlight. 
