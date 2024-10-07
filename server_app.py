# Please excuse the mess. It's prototypical research code. We can clean this up later.
import pandas as pd
import numpy as np
import bokeh
import math
import umap

from bokeh.layouts import column, row
from bokeh.models import Button,Legend,Spacer,HoverTool,ColumnDataSource, Slider, Div, CheckboxGroup,LassoSelectTool, DataTable, TableColumn
from bokeh.plotting import figure, curdoc
import bokeh.palettes as palettes
from bokeh.transform import linear_cmap
from skimage import io
from skimage.util import img_as_ubyte
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

# Genes.csv - List of genes in the feature matrix
# FeatureMatrix.mtx - Feature matrix, stored as a sparse matrix
# ClusterGeneExpression.csv - Gene expression per cell-type cluster
# SpotClusterMembership.csv - Cell type proportions per spot
# SpotPositions.csv - Spot positions and radiuses
# Images/scalefactors_json.json - Scale factors between spot positions
# Images/tissue_hires_image.png - H&E image

spotPostions = pd.read_csv('RedesignChallengeData/SpotPositions.csv')
spotClusterMembership = pd.read_csv('RedesignChallengeData/SpotClusterMembership.csv')
spotPostions['y'] = 17244 - spotPostions['y']
data = pd.merge(spotPostions, spotClusterMembership, on='barcode')
print("loaded data.")

image_path='RedesignChallengeData/Images/tissue_hires_image.png'
he_image_original = io.imread(image_path)
scale_factor=0.11598237
data['LargestCellType'] = spotClusterMembership.iloc[:,4:].idxmax(axis=1)
DotNum = 500
min_value = data.iloc[:, 4:9+4].replace(0, np.nan).min().min()
epsilon = min_value / 10000

# Center log ratio transform
def clr(x, epsilon=epsilon):
    x = x.replace(0, epsilon)
    geometric_mean = np.exp(np.log(x).mean()) 
    return np.log(x / geometric_mean)
clr_data = data.iloc[:, 4:9+4].apply(clr, axis=1)
for i in range(clr_data.shape[1]):
    data[f"clr-{i+1}"] = clr_data.iloc[:,i]
 
# do dimensionality reduction on clr data
def perform_umap():
    dr = umap.UMAP(n_neighbors=30)
    embedding = dr.fit_transform(clr_data.to_numpy())
    data['dr-x'] = embedding[:,0]
    data['dr-y'] = embedding[:,1]

def perform_pca():
    dr = PCA(n_components=2)
    embedding = dr.fit_transform(clr_data.to_numpy())
    data['dr-x'] = embedding[:, 0]
    data['dr-y'] = embedding[:, 1]

perform_pca()

num_clusters = 3
selected_HEclusters = list(range(num_clusters))
kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(clr_data)
data['kmeansCluster'] = kmeans.labels_
centroids = kmeans.cluster_centers_
cellTypes = clr_data.columns.tolist()
palette = palettes.Colorblind[8]*2
# data['color'] = [palette[i] for i in data['kmeansCluster']]
data['clusterName'] = [f'Cluster {i}' for i in data['kmeansCluster']]

kmeansSource = ColumnDataSource(data)
dot_step = data.shape[0]// DotNum
init_dot_idx = data.sort_values(by='dr-x').iloc[[(i + 1) * dot_step - 1 for i in range(DotNum)]].index
DotData = data.loc[init_dot_idx].copy()
DotData['index'] = range(DotNum)
DotDataSource = ColumnDataSource(DotData)
hePlot_visible = True
scatterPlot_visible = True

def HEImage(plot,muted=False):
    global data,num_clusters, imageSource,he_image_original,scale_factor
    he_image = he_image_original.copy()
    alpha = 255
    if muted:
        alpha = 30
    if he_image.shape[2] == 3:  
        alpha_channel = alpha * np.ones((he_image.shape[0], he_image.shape[1], 1), dtype=np.uint8)
        he_image = np.concatenate((he_image, alpha_channel), axis=2)            
    he_image = img_as_ubyte(he_image)
    he_image_flattened = np.flipud(he_image.view(dtype=np.uint32).reshape((he_image.shape[0], he_image.shape[1])))
    dw = he_image.shape[1] / scale_factor
    dh = he_image.shape[0] / scale_factor
    imageSource = ColumnDataSource(data={'image': [he_image_flattened]})
    renderer = plot.image_rgba(image='image', x=0, y=0, dw=dw, dh=dh, source=imageSource)
    return renderer

updateDotBar = True
temp_dots_selected = init_dot_idx

def DotDistributionCallback(attr, old, new):
    global kmeansSource, data, DotDataSource,DotNum, updateDotBar,temp_dots_selected
    global p6
    if not updateDotBar:
        updateDotBar = True
        return  
    DotDataIdx = kmeansSource.selected.indices
    temp_dots_selected = DotDataIdx
    DotData = data.iloc[DotDataIdx].copy()   
    dotNum = DotNum
    if len(DotData) > dotNum:
        part_size = len(DotData) // dotNum
        values = [round(i * part_size)-1 for i in range(1, dotNum + 1)]
        DotData = DotData.sort_values(by='dr-x')
        DotData = DotData.iloc[values]
        temp_dots_selected = DotData.index.tolist()
        temp_dots_selected.sort()
    elif len(DotData) < 1:
        DotData = data.loc[init_dot_idx].copy()
        temp_dots_selected = init_dot_idx
    dotNum = len(DotData)
    DotData.loc[:, 'index'] = range(dotNum)
    DotDataSource.data.update(ColumnDataSource(DotData).data)
    DotDataIdx.sort()

# if a cluster is selected in the checkbox group, update the selected indices of the kmeansSource
def update_selected_clusters(attr, old, new):
    global data, kmeansSource,selected_HEclusters,checkbox_group
    selected_HEclusters = [i for i in checkbox_group.active]
    filtered_data = data[data['kmeansCluster'].isin(selected_HEclusters)]
    selected_indices = filtered_data.index.tolist()
    kmeansSource.selected.indices = selected_indices

def DotSelectedCallback(attr, old, new):
    global data, DotDataSource, DotData, kmeansSource, updateDotBar,temp_dots_selected
    updateDotBar = False
    selected_dot_indices = DotDataSource.selected.indices
    selected_indices = [temp_dots_selected[i] for i in selected_dot_indices]
    kmeansSource.selected.update(indices=selected_indices)

def HESelectionCallback(attr, old, new):
    global data, kmeansSource, imageSource, he_image_original, scale_factor
    he_image = he_image_original.copy()
    SelectedIdx = kmeansSource.selected.indices
    if (hePlot_visible and scatterPlot_visible ) or len(SelectedIdx) == 0:
        alpha_channel = 255 * np.ones((he_image.shape[0], he_image.shape[1], 1), dtype=np.uint8)
        he_image = np.concatenate((he_image, alpha_channel), axis=2)
        he_image = img_as_ubyte(he_image)
        he_image_flattened = np.flipud(he_image.view(dtype=np.uint32).reshape((he_image.shape[0], he_image.shape[1])))
        imageSource.data.update({'image': [he_image_flattened]})
        return
    alpha_channel = 100 * np.ones((he_image.shape[0], he_image.shape[1], 1), dtype=np.uint8)
    he_image = np.concatenate((he_image, alpha_channel), axis=2)
    radius = 142.38 * scale_factor
    radius_squared = radius ** 2
    for index, row in data.iloc[SelectedIdx].iterrows():
        x_center, y_center = row['x'], row['y']
        y_center = 17244 - y_center
        x_center_scaled = x_center * scale_factor
        y_center_scaled = y_center * scale_factor

        x_range = np.arange(int(x_center_scaled - radius), int(x_center_scaled + radius) + 1)
        y_range = np.arange(int(y_center_scaled - radius), int(y_center_scaled + radius) + 1)

        x_grid, y_grid = np.meshgrid(x_range, y_range)
        distance_squared = (x_grid - x_center_scaled) ** 2 + (y_grid - y_center_scaled) ** 2

        mask = distance_squared < radius_squared
        he_image[y_grid[mask], x_grid[mask], 3] = 255

    he_image = img_as_ubyte(he_image)
    he_image_flattened = np.flipud(he_image.view(dtype=np.uint32).reshape((he_image.shape[0], he_image.shape[1])))
    imageSource.data.update({'image': [he_image_flattened]})

def update_clusters(attr, old, new):
    global num_clusters, kmeansSource,data,checkbox_group,selected_HEclusters
    num_clusters = cluster_slider.value  
    # Run k-means clustering with the new number of clusters
    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(clr_data)
    kmeansSource.data['kmeansCluster'] = kmeans.labels_
    kmeansSource.data['clusterName'] = [f'Cluster {i}' for i in kmeans.labels_]
    data['kmeansCluster'] = kmeans.labels_
    new_cluster_options = [f'Cluster {i}' for i in range(num_clusters)]
    checkbox_group.labels = new_cluster_options
    checkbox_group.active = [] 
    selected_HEclusters = list(range(num_clusters))
    # update heimage
    he_image = io.imread(image_path)
    alpha_channel = 255 * np.ones((he_image.shape[0], he_image.shape[1], 1), dtype=np.uint8)
    he_image = np.concatenate((he_image, alpha_channel), axis=2)
    he_image = img_as_ubyte(he_image)
    he_image_flattened = np.flipud(he_image.view(dtype=np.uint32).reshape((he_image.shape[0], he_image.shape[1])))
    imageSource.data.update({'image': [he_image_flattened]})

def store_checkbox_changes(attr, old, new):
    global temp_active
    temp_active = new

def confirm_changes():
    global temp_active
    update_selected_clusters('active', None, temp_active)

lasso_select = LassoSelectTool(continuous = False)
hover = HoverTool(
    tooltips = [
    ("Barcode", "@barcode"),
    ("Cluster", "@kmeansCluster"),
    ("X1", "@X1"),
    ("X2", "@X2"),
    ("X3", "@X3"),
    ("X4", "@X4"),
    ("X5", "@X5"),
    ("X6", "@X6"),
    ("X7", "@X7"),
    ("X8", "@X8"),
    ("X9", "@X9")
    ],
    point_policy="snap_to_data"
)

def configure_circle(plot,fill_alpha=1, line_alpha=1,radius=55):
    renderer = plot.circle('x', 'y', source=kmeansSource,color=mapper, legend_field='clusterName',fill_alpha=fill_alpha, line_alpha=line_alpha, radius=radius)
    plot.add_tools(lasso_select)
    plot.legend.label_text_font_size = '8pt'
    plot.legend.glyph_height = 10
    plot.legend.glyph_width = 10
    plot.legend.spacing = 1
    plot.legend.padding = 1
    return renderer

mapper = linear_cmap('kmeansCluster', palette, 0, len(palette))
tools = "box_select,reset,tap,save,pan,wheel_zoom"
w, h = 600, 550  #415, 415

p1 = figure(width=w, height=h, tools=tools, title='H&E Image with Spots', active_drag="box_select")
heRender = HEImage(p1)
invisible_renderer = p1.circle('x', 'y', source=kmeansSource,color="#00000000", legend_field='clusterName',fill_alpha=0, line_alpha=0, radius=100)
cirRender = configure_circle(p1)
p1.add_tools(HoverTool(
    tooltips = [
    ("Barcode", "@barcode"),
    ("Cluster", "@kmeansCluster"),
    ("X1", "@X1"),
    ("X2", "@X2"),
    ("X3", "@X3"),
    ("X4", "@X4"),
    ("X5", "@X5"),
    ("X6", "@X6"),
    ("X7", "@X7"),
    ("X8", "@X8"),
    ("X9", "@X9"),
    ],
    point_policy="snap_to_data",
    renderers=[cirRender, invisible_renderer]
))

kmeansSource.selected.on_change('indices', DotDistributionCallback, HESelectionCallback)
kmeans.set_selector = lasso_select

p2 = figure(width=w, height=h, title='Similarity of Cell Type Proportions per Spot',
        tools=tools, active_drag="box_select")
#palette = palettes.Set2[num_clusters]
p2.scatter('dr-x', 'dr-y', source=kmeansSource, color=mapper, legend_field='clusterName',
                    muted_color='color', muted_alpha=0.2, size=3)
p2.add_tools(lasso_select, hover)
p2.legend.label_text_font_size = '8pt'
p2.legend.glyph_height = 10
p2.legend.glyph_width = 10
p2.legend.spacing = 1
p2.legend.padding = 1


p3 = figure(width=1550, height=300,title="Cell Type Percentages of Selected Spots",
            toolbar_location='right', tools="box_select,reset,tap,save,xpan,xwheel_zoom", active_drag="box_select")
p3.add_tools(hover)
cellNum = len(cellTypes)
#palette = list(palettes.Light[9])
#palette = list(colorcet.CET_L20[::(256//9)])
palette = []
palette.extend(list(palettes.RdPu[5])[:3])
palette.extend(list(palettes.Blues[4])[:3])
palette.extend(list(palettes.Greens[4])[:3])
palette = palette * math.ceil(cellNum/9)  # Category20c[20] * math.ceil(cellNum/20)
vbars = p3.vbar_stack(cellTypes,x = 'index', width=0.5, color=palette[:cellNum], source=DotDataSource)
p3.x_range.range_padding = 0.02
legend = Legend(items=[(x, [vbars[i]]) for i, x in enumerate(cellTypes)], location="center")
p3.add_layout(legend, 'right')  
p3.legend.title = 'Cell Type' 
DotDataSource.selected.on_change('indices', DotSelectedCallback)

cluster_slider = Slider(start=2, end=8, value=num_clusters, step=1, title="Select number of clusters (k-means)")
cluster_slider.on_change('value_throttled', update_clusters)

temp_active = []

confirm_button = Button(label="Select checked clusters", button_type="success")
cluster_options = [f'Cluster {i}' for i in range(num_clusters)]
selected_clusters_indices = [] 
checkbox_group = CheckboxGroup(labels=cluster_options, active=selected_clusters_indices)

blank = Div(text="  ", width=400, height=15)
checkbox_title = Div(text="Check desired clusters for spot selection", width=400, height=20)
checkbox_group.on_change('active', store_checkbox_changes)
confirm_button.on_click(confirm_changes)

dataTableColumns = [
    TableColumn(field="barcode", title="barcode"),
    TableColumn(field="x", title="x"),
    TableColumn(field="y", title="y"),
    TableColumn(field="X1", title="X1"),
    TableColumn(field="X2", title="X2"),
    TableColumn(field="X3", title="X3"),
    TableColumn(field="X4", title="X4"),
    TableColumn(field="X5", title="X5"),
    TableColumn(field="X6", title="X6"),
    TableColumn(field="X7", title="X7"),
    TableColumn(field="X8", title="X8"),
    TableColumn(field="X9", title="X9"),
    TableColumn(field="kmeansCluster", title="'kmeansCluster")
]
data_table = DataTable(source=kmeansSource, columns=dataTableColumns, width=1000, height=1000)

def toggle_hePlot():
    global hePlot_visible
    hePlot_visible = not hePlot_visible
    update_layout()

def toggle_scatterPlot():
    global scatterPlot_visible
    scatterPlot_visible = not scatterPlot_visible
    update_layout()

def update_layout():
    global hePlot_visible, scatterPlot_visible, p1
    if hePlot_visible and scatterPlot_visible:
        HESelectionCallback('active', None, None)
        heRender.visible = True
        cirRender.visible = True
    elif hePlot_visible:
        HESelectionCallback('active', None, None)
        cirRender.visible = False
    elif scatterPlot_visible:
        heRender.visible = False
    else: 
        hePlot_visible = True
        scatterPlot_visible = True
        HESelectionCallback('active', None, None)
        heRender.visible = True
        cirRender.visible = True

he_button = Button(label="Show/Hide H&E Image", button_type="success")
scatter_button = Button(label="Show/Hide Spot Markers", button_type="success")
he_button.on_click(toggle_hePlot)
scatter_button.on_click(toggle_scatterPlot)

spacer_small = Spacer(width=60, height=40)

# svg rendering (only to create nice looking graphics, slow performance)
#p1.output_backend = "svg"
#p2.output_backend = "svg"
#p3.output_backend = "svg"

# clear_button.width = 200
layout = column(blank,
                cluster_slider,
                checkbox_title,
                row(checkbox_group,spacer_small,confirm_button), 
                spacer_small,
                Div(text="H&E Image plot controls", width=400, height=20),
                he_button, scatter_button
                )
# spacecol = Spacer(width=20, height=415)
tot_layout = row(
                p1,
                Spacer(width=20, height=415),
                p2,
                Spacer(width=20, height=40),
                layout)

curdoc().title = "Visual Compositional Data Analytics for Spatial Transcriptomics"
curdoc().add_root(Div(text="<b>Visual Compositional Data Analytics for Spatial Transcriptomics - David Haegele, Yuxuan Tang, Daniel Weiskopf</b>", width=800, height=20))
curdoc().add_root(tot_layout)
curdoc().add_root(p3)
curdoc().add_root(data_table)
