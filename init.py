# pip install umap-learn
print('Loading in libraries...', end=' ')
import re
import numpy as np
import pandas as pd
import xgboost as xgb
import seaborn as sns
import umap.umap_ as umap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.stats import uniform
from sklearn.svm import SVR, NuSVR
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.linear_model import LinearRegression, RidgeCV, LassoCV, Lasso, Ridge, ElasticNet, ElasticNetCV
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, train_test_split

print('done.', end='\n')

print('Loading in functions...', end=' ')
def standardise(Dataframe):
    # standardises all the columns of a Dataframe
    standardised = (Dataframe - Dataframe.mean())/Dataframe.std()
    standardised =  standardised.fillna(0)
    return(standardised)

def validation_plot(y_val, y_pred, y_axis_label='age'):
    # if y_val is a dataframe:
    if isinstance(y_val, pd.DataFrame):
        y_val = y_val.values

    r2 = r2_score(y_val, y_pred)
    diff = np.round(np.mean(np.abs(y_val-y_pred)), 2)
    # sort data because it's easier to interpret
    zipped = zip(y_pred, y_val)
    sorted_zipped = sorted(zipped, key=lambda x: x[1])
    y_pred, y_val = zip(*sorted_zipped)
    # plot
    x = range(0,np.size(y_pred))
    fig, ax = plt.subplots(figsize=(17, 6), dpi=400)
    ax.set_title(f'R-squared: {r2}\nMean Absolute Difference: {diff}')
    ax.set_xlabel('participant');
    ax.set_ylabel(y_axis_label);
    ax.scatter(x, y_val);
    ax.scatter(x, y_pred, color='orange');
    ax.legend(['Observation', 'Predicted']);

def correlation_list(predictor, response, analyte_conv=None):
    colnames = predictor.columns
    p_vec = []

    for i in range(len(colnames)):
        pearson_corr = np.corrcoef(response, predictor.iloc[:, i])[0, 1]
        p_vec.append(pearson_corr)

    df = pd.DataFrame({'Pearson Correlation':p_vec, 'SeqIndex':colnames})
    df = df.iloc[np.abs(df['Pearson Correlation']).argsort()[::-1]]

    if analyte_conv is not None:
        df = df.merge(analyte_conv)

    return df

def search_metadata(metadata, search_term):
    # metadata = R2_Meta
    # search_term = r'R1|R2'

    # Perform search
    matches = [col for col in metadata.columns if re.search(search_term, col, re.IGNORECASE)]

    matching = metadata[matches]

    # Print results
    print(f'Found {len(matches)} matches')
    print('Match   |   contains nan values')
    print(metadata[matches].isna().sum())

    out = matching.head(3).T
    # set column names
    out.columns = ['Sample 1', 'Sample 2', 'Sample 3']
    return out

def minmax_scale(df):
    scaler = MinMaxScaler()
    scaled_data = scaler.fit_transform(df)

    if isinstance(df, pd.DataFrame):
        scaled_data = pd.DataFrame(scaled_data, columns=df.columns)

    return scaled_data

results = pd.DataFrame()    # initialise results dataframe

def hyperparameter_search_EN(dataloader, param_grid, results, clock_name='EN', cores=-1):
    # Unpack dataloader
    x_train = dataloader['x_train']
    y_train = dataloader['y_train']
    x_test = dataloader['x_test']
    y_test = dataloader['y_test']


    # Create Elastic Net regression model
    elastic_net = ElasticNet()

    # Perform hyperparameter tuning using GridSearchCV
    grid_search = GridSearchCV(estimator=elastic_net, param_grid=param_grid, cv=5, scoring='neg_mean_squared_error', n_jobs=cores)
    grid_search.fit(x_train, y_train)

    # Get the best hyperparameters
    best_params = grid_search.best_params_
    print("Best Hyperparameters:", best_params)

    # Use the best model for prediction
    best_model = grid_search.best_estimator_
    y_pred_train = best_model.predict(x_train)
    y_pred_test = best_model.predict(x_test)

    # Evaluate the model
    mse_train = mean_squared_error(y_train, y_pred_train)
    mse_test = mean_squared_error(y_test, y_pred_test)
    r2_train = r2_score(y_train, y_pred_train)
    r2_test = r2_score(y_test, y_pred_test)

    print("Mean Squared Error (Train):", mse_train)
    print("Mean Squared Error (Test):", mse_test)
    print("R-squared (Train):", r2_train)
    print("R-squared (Test):", r2_test)

    # Append results to results dataframe
    output = ({
        'model': clock_name,
        'hyperparameters': [best_params],
        'mse_train': mse_train,
        'mse_test': mse_test,
        'r2_train': r2_train,
        'r2_test': r2_test
    })

    results = pd.concat([results, pd.DataFrame(output)], ignore_index=True)

    return results, best_model, y_pred_test

def load_in_data(dataset):
    R1_path = 'data/'+dataset+'_phase1.csv'
    R2_path = 'data/'+dataset+'_phase2.csv'
    probe_path = 'data/probe_metadata.csv'

    # Round 1
    try:
        print('Loading in round 1 data...')
        R1_SOMAMeta_raw = pd.read_csv(R1_path, low_memory=False)

    except:
        print('Round 1 data not found...')

    # Round 2
    try:
        print('Loading in round 2 data...')
        R2_SOMAMeta = pd.read_csv(R2_path, low_memory=False)
    except:
        print('Round 2 data not found...')

    # TargetID to EntrezID
    try:
        print('Loading in probe metadata...')
        analyte = pd.read_csv(probe_path)
        analyte_conv = analyte[['SeqIndex', 'EntrezGeneSymbol']]
    except:
        print('Probe metadata not found...')


    print('Cleaning data...')

    # Round 1
    R1_SOMAMeta_raw['VisitB_Date'] = pd.to_datetime(R1_SOMAMeta_raw['VisitB_Date'])
    R1_SOMAMeta_raw['year-month'] = R1_SOMAMeta_raw['VisitB_Date'].dt.strftime('%Y-%m')

    # Remove batch effected months April, May, June, July (2003)
    months_to_remove = [4, 5, 6, 7] 
    R1_SOMAMeta = R1_SOMAMeta_raw.copy()
    R1_SOMAMeta = R1_SOMAMeta[~((R1_SOMAMeta['VisitB_Date'].dt.year == 2003) & 
                                (R1_SOMAMeta['VisitB_Date'].dt.month.isin(months_to_remove)))]
    R1_SOMAMeta['year-month'] = R1_SOMAMeta['VisitB_Date'].dt.strftime('%Y-%m')

    # Round 2
    R2_SOMAMeta['Examination date R2'] = pd.to_datetime(R2_SOMAMeta['Examination date R2'])
    R2_SOMAMeta['year-month'] = R2_SOMAMeta['Examination date R2'].dt.strftime('%Y-%m')


    # Set Index
    R1_SOMAMeta = R1_SOMAMeta.set_index(R1_SOMAMeta['visit_id'])
    R2_SOMAMeta = R2_SOMAMeta.set_index(R2_SOMAMeta['visit_id'])

    # Select proteins
    R1_SOMA = R1_SOMAMeta.loc[:, R1_SOMAMeta.columns.str.startswith('seq')]
    R2_SOMA = R2_SOMAMeta.loc[:, R2_SOMAMeta.columns.str.startswith('seq')]

    # there is one protein in R2 that is not in R1 (seq.20367.6)
    shared_proteins = list(set(R1_SOMA.columns) & set(R2_SOMA.columns)) 
    R1_SOMA = R1_SOMA.loc[:, shared_proteins]
    R2_SOMA = R2_SOMA.loc[:, shared_proteins]


    # Select metadata
    R1_Meta = R1_SOMAMeta.loc[:, ~R1_SOMAMeta.columns.str.startswith('seq')]
    R2_Meta = R2_SOMAMeta.loc[:, ~R2_SOMAMeta.columns.str.startswith('seq')]

    # Select age data
    R1_Age = R1_SOMAMeta['Age_years'].to_frame().rename(columns={'Age_years': 'Age'})
    R2_Age = R2_SOMAMeta['Age R2 (year)'].to_frame().rename(columns={'Age R2 (year)': 'Age'})

    print('Done.')
 
    return R1_SOMAMeta_raw, R1_SOMAMeta, R2_SOMAMeta, R1_SOMA, R2_SOMA, R1_Meta, R2_Meta, R1_Age, R2_Age, analyte_conv


print('done.')

