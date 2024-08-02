#!/usr/bin/env python3

import os
from pathlib import Path
import pandas as pd
import dash
from dash import dash_table
from dash import dcc
from dash import html
import plotly.graph_objs as go
import plotly
import plotly.figure_factory as ff
import time
import typer

from functools import partial

from dash.dependencies import Input, Output

from typing import List

from chemdash.chemistry import df_canonicalize_from_smiles, df_smiles_to_rdkit_mol, \
    df_get_image_file_url, df_add_ecfp_1024_4_fps, df_add_maccs_fps

from chemdash.molecular_descriptors import MolecularDescriptors
from chemdash.dimensionality import separate_metadata, generate_lap_grid, \
    add_dummy_rows, generate_tsne, generate_umap

from chemdash.multiprocess import parallelize_dataframe_func


MAX_COMPOUNDS = 10000


def prepare_data(df, fp_type='maccs', plot_type='umap'):
    print('canonicalizing smiles')
    t0 = time.time()
    print(f"Start time: {t0}")
    df['canon_smiles'] = df_canonicalize_from_smiles(df, 'smiles')
    print('converting to mols')
    t1 = time.time(); e1 = t1 - t0
    print(f"elapsed: {e1} ({t1})")
    df['mol'] = df_smiles_to_rdkit_mol(df, 'canon_smiles')
    print('generating images')
    t2 = time.time(); e2 = t2 - t1
    print(f"{e2} ({t2}")
    print(f"elapsed: {e2} ({t2})")

    df['Structure'] = parallelize_dataframe_func(df, partial(
        df_get_image_file_url, mol_col='mol', id_col='Compound_id'))

    if fp_type == 'maccs':
        print ('creating fingerprints')
        print(time.time())
        df, bit_columns, bit_column_prefix = df_add_maccs_fps(df, 'mol')
    else:
        #will run an ecfp 1024 radius 2
        df, bit_columns, bit_column_prefix = df_add_ecfp_1024_4_fps(df, 'mol')

    #create Molecular Descriptor object and calc descriptors on the df
    print ('calc mol descriptors')
    print(time.time())
    t3 = time.time(); e3 = t3 - t2
    print(f"elapsed: {e3} ({t3})")

    md = MolecularDescriptors(mol_field='mol', smiles_field='canon_smiles', id_field='Compound_id')
    df = md.generate_descriptors(df)

    #remove the mol column/smiles column from import--keeping the canon smiles
    df = df.drop(columns=['mol', 'smiles'], axis=1)

    #prepare for generating plots
    df = add_dummy_rows(df, 'Compound_id', bit_columns)
    x, metadata = separate_metadata(df, bit_column_prefix)

    if plot_type == 'tsne':
        df = generate_tsne(df, x, metadata)
        df = generate_lap_grid(df, 'tsne_x', 'tsne_y')
    else:
        print ('creating umap')
        print(time.time())
        df = generate_umap(df, x, metadata)
        print ('create lap grid')
        print(time.time())
        df = generate_lap_grid(df, 'umap_x', 'umap_y')

    df = df[df['Compound_id'] != 'Placeholder'].reset_index(drop=True)

    print('Done Cleaning')
    t4 = time.time(); e4 = t4 - t3
    print(f"{e4} ({t4}")
    print(f"elapsed: {e4} ({t4})")

    return df


def get_table_df(df:pd.DataFrame)->(pd.DataFrame, List):
    cols_to_remove = ['tsne_x', 'tsne_y', 'lap_x', 'lap_y', 'umap_x', 'umap_y',]
    cols_to_keep =  [x for x in df.columns if x not in cols_to_remove]
    df_table = df[cols_to_keep]
    #want to drop additional columns from the list which will be displayed
    cols_to_keep.remove('Structure')
    cols_to_keep.remove('canon_smiles')
    cols_to_keep.remove('Compound_id')
    cols_to_keep.sort()

    return df_table, cols_to_keep


def get_table_layout(df):
    df_table, df_table_cols = get_table_df(df)

    table_layout = dash_table.DataTable(
        id='data_table',
        columns=[
            {   #want the structue col first
                'id': 'Structure',
                'name': 'Structure',
                'presentation': 'markdown',
            },
            {   'id': 'Compound_id',
                'name': 'Compound_id',
            },
            {   'id': 'canon_smiles',
                'name': 'canon_smiles',
            },] +

            [{"name": i, "id": i, "selectable": True} for i in df_table_cols]
            ,
        #unfortuantely there are some serious formatting issues when the columns get fixed...
        #fixed_columns={ 'headers': True, 'data': 1},
        data=df.to_dict('records'),
        #filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        #row_selectable="multi", #turned this off until figure out virtual selected rows
        selected_columns=['MolecularWeight'],
        #selected_rows=[],
        page_current=0,
        page_size=25,
        page_action='native',
        style_cell={'width': '150px'},
        #editable=True,
        style_table={'overflowX': 'scroll',
                     'maxHeight': '700px',
                     'overflowY': 'scroll'
                     },
        fixed_rows={'headers': True, 'data': 0},
        style_data_conditional=[
        {
            'if': {'row_index': 'odd'},
            'backgroundColor': 'rgb(248, 248, 248)'
        }],
    )

    return table_layout


def get_sc_figure(df, colorby, x, y):
    figure = {
        'data': [
            go.Scattergl(
                x=df[x],
                y=df[y],
                mode='markers',
                opacity=0.9,
                marker={

                    'color': df[colorby],
                    'line': {'width': 0.5, 'color': 'white'},
                    'symbol': 'square',
                    'colorbar': dict(thickness=20)
                },
                name="sarm_data",
                marker_colorscale=plotly.colors.sequential.RdBu_r,

            ),
        ],
        'layout': go.Layout(
            height=600,
            width=800,
            xaxis={
                   'showline': False,
                   'showgrid': False,
                   'showticklabels': False},
            yaxis={
                   'showline': False,
                   'showgrid': False,
                   'showticklabels': False},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 1, 'y': 1},
            hovermode=False,
            dragmode='select',
            plot_bgcolor='Black',

        )
    }
    return figure

def get_mgm_graph():

        mgm_graph = dcc.Graph(
            id='mgm_graph',
            #config={'displayModeBar': True,},
        )

        return mgm_graph


def get_tsne_graph():
    tsne_graph = dcc.Graph(
        id='tsne_graph',
        #config={'displayModeBar': True, },
    )

    return tsne_graph


def get_umap_graph():
    umap_graph = dcc.Graph(
        id='umap_graph',
        #config={'displayModeBar': True, },
    )

    return umap_graph


def get_columns_dropdown(df: pd.DataFrame):
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    all_cols = df.select_dtypes(include=numerics).columns
    cols_to_remove = ['tsne_x', 'tsne_y', 'lap_x', 'lap_y', 'id', 'umap_x', 'umap_y',]
    cols_to_keep = [x for x in all_cols if x not in cols_to_remove]
    cols_to_keep.sort()
    return cols_to_keep


def create_chemdash_app(df: pd.DataFrame):
    app = dash.Dash(__name__, assets_folder='assets')

    #useful when trying to get callbacks sorted out...
    #app.config.suppress_callback_exceptions = True

    #layout...
    app.layout = html.Div([html.A([
                                html.Img(src='assets/CDlogo_color.png',
                                            style = {
                                            'height': '10%',
                                            'width': '10%',
                                            'float': 'right',
                                            'position': 'relative',
                                            'padding-top': 0,
                                            'padding-right': 0
                                            })
                                        ],
                                    href='https://www.cognitivedataworks.com'),
        html.H4("GAN Compound Explorer"),
        html.Div([
        html.Div([
            html.Div(
                [dcc.Dropdown(
                        id='color_by_dropdown',
                        options=[
                            {"label": i, "value": i} for i in get_columns_dropdown(df)
                            ],
                        value='MolecularWeight', clearable=False
                        ),
                dcc.Tabs(id='tabs', children=[
                dcc.Tab(label='MGM', id='mgm_tab', children=[get_mgm_graph()],
                        ),
                dcc.Tab(label='UMAP', id='umap_tab', children=[get_umap_graph()],
                        ),
                ]),
            ])],  className="six columns"),
        html.Div([dcc.Graph(id='density-plot')], className="six columns")], className="row"),
        html.Div([get_table_layout(df),
                 ])
    ])





    @app.callback(
        [Output('mgm_graph', 'figure'),
         Output('umap_graph', 'figure')],
        [Input('color_by_dropdown', 'value')]
    )
    def generate_scatter_figs(colorby):
        return get_sc_figure(df, colorby, 'lap_x', 'lap_y'), get_sc_figure(df, colorby, 'umap_x', 'umap_y')

    #bit of hack in order to let you switch between plots and select data
    #without first clearing the data on previous plot
    @app.callback(
        [Output('mgm_graph','selectedData'),
         Output('umap_graph','selectedData')],
        [Input('tabs', 'value'),
        Input('mgm_graph', 'clickData'),
        Input('umap_graph', 'clickData')
    ]
    )
    def update_tab_select(value, click1, click2):
        print (value, click1, click2)
        return {}, {}

    @app.callback(
        Output('data_table', 'data'),
        [Input('mgm_graph', 'selectedData'),
        Input('umap_graph', 'selectedData'),
    ]
    )
    def update_table(selectedDataSarm, selectedDataUmap):

        if selectedDataSarm:
            if len(selectedDataSarm['points']) == 0 :
                return df.to_dict('records')
            match_idx = [x['pointIndex'] for x in selectedDataSarm['points']]
            return df.iloc[match_idx].to_dict('records')
        elif selectedDataUmap:
            if len(selectedDataUmap['points']) == 0:
                return df.to_dict('records')
            match_idx = [x['pointIndex'] for x in selectedDataUmap['points']]
            return df.iloc[match_idx].to_dict('records')
        else:
            return df.to_dict('records')



    @app.callback(
        [Output('data_table', 'style_data_conditional'),
        Output('density-plot', 'figure')],
        [Input('data_table', 'selected_columns'),
         Input('data_table', 'data'),]
    )
    def get_dist_plot(selected_columns, data):

        col = str(selected_columns[0])
        if df.shape[0] == len(data):
            hist_data = [df[col].dropna()]
            group_labels = ['All Compounds']
            bin_size = max(df[col]) / 100
        else:
            mol_index = [i['id'] for i in data]
            hist_all = df[~df['id'].isin(mol_index)][col].dropna()
            hist_selected = df[df['id'].isin(mol_index)][col].dropna()
            hist_data = [hist_all, hist_selected]
            group_labels = ['All Compounds', 'Currently_Selected']
            bin_size = max(df[col]) / 100

        #known 'bug' here you need to select more than one data point from the SARM/TSNE in order to
        #create the updated distplot.  I did not have time to look into why this is
        #histnorm can be ['', 'percent', 'probability', 'density', 'probability density']
        fig =  ff.create_distplot(hist_data, group_labels, histnorm = 'density', bin_size=bin_size)

        fig.layout.update(title=col, height=600, width=800)

        select_cols = [{
            'if': {'column_id': i},
            'background_color': '#D2F3FF'
        } for i in selected_columns]

        return select_cols, fig

    return app


def main(
    input_file: Path,
    port: int = 8000,
    debug: bool = False
):
    """
    Run the chemdash server.

    :param input_file: Path to the input file of compounds.
    :param port: Port number to use (default is 8000).
    :param debug: Debug mode (default is False).
    """
    df = pd.read_csv(str(input_file))
    if not all(x in df.columns for x in ['smiles', 'Compound_id']):
        raise ValueError('Your dataset needs to contain "smiles" and "Compound_id" columns!')
    if df.shape[0] > MAX_COMPOUNDS:
        raise ValueError('This tool supports a maxiumum of 10k compounds!')

    typer.echo('There are ' + str(df.shape[0]) + ' compounds in your dataset')

    # create asset directory for mol image files
    if not os.path.isdir("assets/mol_images"):
        os.mkdir('assets/mol_images')

    df = prepare_data(df)
    df['id'] = df.index

    app = create_chemdash_app(df)

    typer.echo(f"Launching server: http://0.0.0.0:{port}, debug: {debug}")
    app.run_server(host="0.0.0.0", port=port, debug=debug)
    typer.echo(f"Server exited.")


if __name__ == '__main__':
    # NOTE: THe parallelization logic in prepare() will re-load this module in sub-processes,
    # All "real work" is done in this block, with the main module body doing only definitions!
    typer.run(main)

