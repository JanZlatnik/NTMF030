import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator



# Nastavení základních vlastností grafu
name = '2.cross_1'
aspect_ratio = 'square'
filenames = [r"D:\Univerzita\Quantum scattering theory\Potential_scattering\Task 2\cross_sections.txt"]
x_column = 1
y_columns = [3]
show_legend = True  
legend_title = 0
legend_frame = False 

# Nastavení velikostí
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

# Legendy, popisky, barvy a styly čar s escape sekvencemi pro LaTeX
X_axis = r'Energy$\,(\mathrm{Ha})$'
Y_axis = r'Cross section$\,(a_0^2)$'
legends = [[r"$l=1$"],[0],[0],[0]]
labels = [[0,0,0,0]]
labels_position = [[(2.63,10**(-1.7)),(1.4, 10**(-1.5)),(1.4, 10**(-3.1)),(2.2, 10**(-1.03))]]  
colors = [['green','red','blue','purple']]
markers = []
line_styles = [['-','--','-.',':']]
line_widths = [["",'','','']]
ylog = False
xlog = False

# Nastavení pro omezení os a krokování
limit_x = False
limit_y = True 
x_range = [0, 50]  
y_range = [0, 100] 
minstep_x = False
minstep_y = False 
step_x = 0.5  
step_y = 10 


# Nastavení LaTeX formátování a fontů
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Latin Modern Roman']

plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = ["Latin Modern Roman"]
plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)          
plt.rc('axes', titlesize=MEDIUM_SIZE)    
plt.rc('axes', labelsize=MEDIUM_SIZE)    
plt.rc('xtick', labelsize=SMALL_SIZE)    
plt.rc('ytick', labelsize=SMALL_SIZE)    
plt.rc('legend', fontsize=LEGEND_SIZE)       
plt.rc('figure', titlesize=BIGGER_SIZE) 

# Nastavení velikosti a poměru stran figury
if aspect_ratio == 'golden_ratio':
    golden_ratio = (1 + 5 ** 0.5) / 2
    fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
else:
    fig, ax = plt.subplots(figsize=(SIZE, SIZE))

# Načtení a vykreslení dat z více souborů
for file_idx, filename in enumerate(filenames):
    data = np.loadtxt(filename, comments='#', delimiter=None)
    x = data[:, x_column-1]
    for col_idx, y_col in enumerate(y_columns):
        y = data[:, y_col-1]
        label = legends[file_idx][col_idx]
        y = y
        x_filtered = x
        if label:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x_filtered, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x_filtered, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        else:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x_filtered, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x_filtered, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        if labels[file_idx][col_idx] != 0:
            abs_x_pos = labels_position[file_idx][col_idx][0]
            abs_y_pos = labels_position[file_idx][col_idx][1]
            ax.text(abs_x_pos, abs_y_pos, labels[file_idx][col_idx], color=colors[file_idx][col_idx], verticalalignment='bottom')

# Nastavení vzhledu grafu
ax.set_xlabel(X_axis)
ax.set_ylabel(Y_axis)
if limit_x:
    ax.set_xlim(x_range)
if limit_y:
    ax.set_ylim(y_range)
ax.set_title('')
if show_legend:
    if legend_title != 0:
        ax.legend(frameon=legend_frame, title=legend_title, loc = "upper right")
    else:
        ax.legend(frameon=legend_frame)
ax.tick_params(axis='x', direction='in', which='both', top=True, bottom=True, labelbottom=True, length = 6)
ax.tick_params(axis='y', direction='in', which='both', left=True, right=True, labelleft=True, length = 6)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='x', direction='in', which='minor', top=True, bottom=True, labelbottom=True, length = 2)
ax.tick_params(axis='y', direction='in', which='minor', left=True, right=True, labelleft=True, length = 2)
ax.grid(False)
if ylog:
    ax.set_yscale('log')
if xlog:
    ax.set_xscale('log')
if limit_x:
    ax.set_xlim(x_range)
    if minstep_x:
        ax.xaxis.set_major_locator(plt.MultipleLocator(step_x))
if limit_y:
    ax.set_ylim(y_range)
    if minstep_y:
        ax.yaxis.set_major_locator(plt.MultipleLocator(step_y))

fig.savefig(f'{name}.pdf', format='pdf')
plt.show()
