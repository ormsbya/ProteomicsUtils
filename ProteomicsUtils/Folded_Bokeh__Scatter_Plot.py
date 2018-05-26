from bokeh.plotting import figure, show, output_file
from bokeh.palettes import brewer
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar,
    glyphs,
    Span
)
from bokeh.io import export_svgs
from bokeh.layouts import gridplot
import pandas as pd
import os
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info("Import successful")


def bokeh_scatter_maker(df, x_col, y_col, c_col, title, hover_list, to_svg=False):
    """Generates a Bokeh figure object for the supplied data, with hover interactivity. Points are coloured according to a datacolumn
    Parameters
    ----------
    df : dataframe
        Pandas DataFrame containing the data to be plotted (x, y, and colours, plus optional hover information) as columns.
    x_col : str
        Column name for x_data, also used to label the x-axis
    y_col : str
        Column name for y_data, also used to label the y-axis
    c_col : str
        Column name for data used to set colour of points
    title : str
        Title of figure to be generated
    hover_list : list of tuples
        Each tuples contains the name of the parameter to be shown on hover (e.g. Protein), with the column to be used for mapping (e.g. @{Master Protein Accessions}), both as strings. Column names with spaces should be encased in {}, and map columns preceded with @.
    to_svg : bool, False (default)
        If true, the save button on the generated html plot will save an svg. Be aware this can affect the interactive functions in html version.
    Returns
    -------
    fig
        Bokeh figure object, which can be plotted (using show(fig)) or added to grid layout.
    """

    source = ColumnDataSource(df)
    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

    colors = list(reversed(brewer['Reds'][9]))#brewer['RdYlBu'][25]#["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=colors, low=df[c_col].min(), high=df[c_col].max())

    fig = figure(title=title, plot_width=500, plot_height=500, tools=TOOLS, toolbar_location='below')

    #creating objects to be added to the figure
    vline = Span(location=0, dimension='height', line_color='red', line_width=1, line_dash='dotted')
    hline = Span(location=0, dimension='width', line_color='red', line_width=1, line_dash='dotted')
    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                     ticker=BasicTicker(desired_num_ticks=len(colors)),
                     label_standoff=6, border_line_color=None, location=(0, 0))

    #adding all elements to the figure plot
    fig.renderers.extend([vline, hline])
    fig.grid.grid_line_color = None
    fig.background_fill_color = None
    fig.xaxis.axis_label = str(x_col)
    fig.yaxis.axis_label = str(y_col)
    fig.scatter(x=x_col,
          y=y_col,
          marker='circle', size=15,
          source=source,
          fill_color={'field': c_col, 'transform': mapper},
              line_color="navy", alpha=0.5)
    fig.add_layout(color_bar, 'right')
    fig.select_one(HoverTool).tooltips = hover_list
    #fig.select_one(HoverTool).formatters={'Gene name' : 'printf', 'Ontology' : 'printf',# use 'printf' formatter}
    if to_svg:
        fig.output_backend = "svg"

    return fig

def main(input, output_path, sample_name):

    #Checking input type: expects df or file path
    if isinstance(input, pd.DataFrame):
        logger.info(f"Input detected as DataFrame.")
        df = input
    elif os.path.isfile(input):
        logger.info(f"Collecting data from {input}...")
        df = pd.read_excel(input)
    else:
        logger.info(f"Incorrect input format detected. Please pass full file path, or a dataframe as input.")

    logger.info(f"Preview data: {df.head(10)}")

    output_file(output_path, title=sample_name+"_scatter")
    logger.info(f"Output html will be saved to {output_path}")

    hovers = [('Protein', '@{Master Protein Accessions}'),
        ('Gene Name', '@{Gene names}'),
        ('Ontology', '@Ontology'),]

    fig = bokeh_scatter_maker(df=df, c_col='-Log10 p-Value', y_col='Log2 Average Ratio', x_col='Log2 Average Non-cys', title=sample_name, hover_list=hovers)
    show(fig)

    return fig

if __name__ == "__main__":
    ## Default data for main function
    input_path='C:/Users/dezer_000/Desktop/dezeraecox.com/Plotting Post - TPE-Tunicamycin/Data_to_plot_gene_info_mapped.xlsx'
    output_path = "C:/Users/dezer_000/Desktop/dezeraecox.com/Plotting Post - TPE-Tunicamycin/grid_scatter.html"
    sample_name='Tunicamycin'
    main(input_path, output_path, sample_name)
