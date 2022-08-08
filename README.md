# FlaGs Results Visualisation

Complex programs can fail to deliver their results without an adequate visualization interface. FlaGs (Flanking Genes) is a bioinformatics tool that detects homologous gene clusters and outputs a graphical visualization along with a phylogenetic tree. Recent updates to the FlaGs algorithm has made this tool domain-aware, however it is lacking the visual interface to support such feature. 

This project aims to:
1.	improve the original plotting script to avoid some previous limitations of the program (using matplotlib, matplotlib.py), 
2.	develop a semi-interactive, easily sharable FlaGs result suitable for a limited number of domains per protein (using plotly, semi-interactive_html.py), and
3.	develop a fully interactive FlaGs application in order to provide complete and convenient access to all domain data (using Plotly Dash, interactive_dash.py).



## Contents

Overview of repository files including the script and each files it requires in a share directory. All files were acquired using FlaGs2, https://github.com/GCA-VH-lab/FlaGs2 
![image](https://user-images.githubusercontent.com/100831180/183397550-e7fa4cac-fa95-4514-b5c2-21b08d68ce36.png)
	

## How to Install and Run the Project

Downloading the script with its corresponding file/s in a shared directory should be enough to run and display the results locally. Some dependencies, such as dash_bio might need to be downloaded in order to run the interactive_dash.py script. 

The matplotlib.py will display the results using matplotlib viewing window once the code is run. The semi-interactive_html.py script will open a HTML server to display results. The interactive_dash.py script requires a server to run, running it locally requires opening the link pasted at the top of the script to open the results in an interactive dashboard. 


## How to Use the Project
1.	The matplotlib.py script bypasses some limitations the FlaGs1 visualization had. This script is already incorporated in the follow-up version, FlaGs2. 
2.	The semi-interactive_html.py script allows the user to visualize domains inside the genes as discovered by the domains search option of FlaGs2. It allows limited interactivity whereby the user can toggle domains and proteins using the legend on the right-hand side. Given the simplicity and limitations to this layout, it is recommended to use a limited number of domains per proteins (e.g. 3-4 domains). This is because using more makes the legend less convenient for use due to its length. 
3.	The full interactive_dash.py script allows full interactivity and access to all domains regardless of the number of domains. It includes a domain plot which is created upon the selection of a gene as well as a sequence alignment of the query proteins sequences (black arrows) aligned with the MAFFT multiple sequence alignment method.

## Credits
This project was made to improve the bioinformatic tool FlaGs developed by the Atkinson-Hauryliuk lab, https://github.com/GCA-VH-lab.



