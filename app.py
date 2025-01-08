import dash
from dash import dcc, html, Input, Output, State, dash_table
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import time  

df = pd.read_csv('/Users/u1536635/Downloads/notebooks/high_qed_coconut_desc.csv')  # Ensure this file includes the necessary descriptors

numeric_columns = [
    'MolecularWeight', 'LogP', 'TPSA', 'NumRotatableBonds', 'NumAromaticRings',
    'HBondDonors', 'HBondAcceptors', 'HeavyAtomCount', 'MolFractionCSP3',
    'RingCount', 'Chi0', 'Chi1', 'NPLikeness', 'QED_Drug_Likeness', 'SyntheticAccessibilityScore', 'Cluster'
]
df[numeric_columns] = df[numeric_columns].round(2)
df = df[numeric_columns + ['COCONUT_ID', 'SMILES']].dropna()

app = dash.Dash(__name__, suppress_callback_exceptions=True)
# app.title = "COCONUT Similarity Search"

app.layout = html.Div([
    html.H1("COCONUT Database", style={'text-align': 'center'}),
            html.H1("Tanimoto Similarity Search", style={'textAlign': 'center'}),
            html.Div([
                html.Label("Enter a SMILES string:"),
                dcc.Input(id='smiles-input', type='text', style={'width': '50%'}, placeholder="Enter SMILES"),
                html.Button("Search", id='search-button', n_clicks=0)
            ], style={'textAlign': 'center', 'marginBottom': '20px'}),
            html.Div([
                html.Label("Time taken for search:"),
                html.Div(id='time-taken', style={'fontWeight': 'bold'})
            ], style={'textAlign': 'center', 'marginBottom': '20px'}),
            dash_table.DataTable(
                id='similarity-table',
                columns=[  
                    {'name': 'COCONUT_ID', 'id': 'COCONUT_ID'},
                    {'name': 'SMILES', 'id': 'SMILES'},
                    {'name': 'MolecularWeight', 'id': 'MolecularWeight'},
                    {'name': 'LogP', 'id': 'LogP'},
                    {'name': 'Cluster', 'id': 'Cluster'},
                    {'name': 'SimilarityScore', 'id': 'SimilarityScore'}
                ],
                style_table={'height': '300px', 'overflowY': 'auto'},
            )
        ])

# Callback for Similarity Search
@app.callback(
    [Output('similarity-table', 'data'),
     Output('time-taken', 'children')],
    Input('search-button', 'n_clicks'),
    State('smiles-input', 'value')
)
def perform_similarity_search(n_clicks, input_smiles):
    if not input_smiles:
        return [], "N/A"

    try:
        start_time = time.time()
        query_mol = Chem.MolFromSmiles(input_smiles)
        if query_mol is None:
            return [{"COCONUT_ID": "Invalid SMILES", "SMILES": "N/A", "MolecularWeight": "N/A", "LogP": "N/A", "Cluster": "N/A", "SimilarityScore": "N/A"}], "N/A"

        query_fp = FingerprintMols.FingerprintMol(query_mol)

        # Calculate similarity for each molecule in the dataset
        similarities = []
        for idx, row in df.iterrows():
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol:
                mol_fp = FingerprintMols.FingerprintMol(mol)
                similarity = DataStructs.FingerprintSimilarity(query_fp, mol_fp)
                similarities.append((similarity, row['COCONUT_ID'], row['SMILES'], row['MolecularWeight'], row['LogP'], row['Cluster']))

        similarities = sorted(similarities, key=lambda x: x[0], reverse=True)
        top_similar = [{
            "COCONUT_ID": sim[1],
            "SMILES": sim[2],
            "MolecularWeight": sim[3],
            "LogP": sim[4],
            "Cluster": sim[5],
            "SimilarityScore": f"{sim[0]:.4f}"  
        } for sim in similarities[:10]]
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        time_taken = f"{elapsed_time:.2f} minutes"

        return top_similar, time_taken

    except Exception as e:
        return [{"COCONUT_ID": "Error", "SMILES": str(e), "MolecularWeight": "N/A", "LogP": "N/A", "Cluster": "N/A", "SimilarityScore": "N/A"}], "N/A"

if __name__ == '__main__':
    app.run_server(debug=True)
