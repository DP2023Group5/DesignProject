from scipy.stats import pearsonr
from statsmodels.stats import multitest
import statsmodels.api as sm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def process_outputs(y_test, y_pred):
    y_pred = pd.DataFrame(y_pred, columns=['Predicted Age'])
    y_pred.index = y_test.index

    test_corr_data = pd.concat([y_test, y_pred], axis=1)

    test_corr_data['Error'] = test_corr_data['Predicted Age'] - test_corr_data['Age']

    print('MAE:', np.mean(abs(test_corr_data['Error'])))

    return test_corr_data

def fit_smoother(test_corr_data, show_plot=True):

    # Fit a local regression (LOESS)
    lowess = sm.nonparametric.lowess
    frac = 2/3  # Fraction parameter set to 2/3

    # Sorting data by chronological age for smoother plot
    test_corr_data_sorted = test_corr_data.sort_values(by='Age')

    # Fit the lowess model
    predicted_smooth = lowess(test_corr_data_sorted['Predicted Age'], 
                            test_corr_data_sorted['Age'], 
                            frac=frac)

    # Extracting smoothed values
    smoothed_values = predicted_smooth[:, 1]

    individual_age_gaps = test_corr_data_sorted['Predicted Age'] - smoothed_values

    # Reindex individual_age_gaps to match the original index of test_corr_data
    individual_age_gaps = individual_age_gaps.reindex(test_corr_data.index)

    # Add individual age gaps to the DataFrame
    test_corr_data['Individual age gaps'] = individual_age_gaps

    if show_plot:
        # Plotting true age vs error
        plt.figure(figsize=(10, 6))
        plt.scatter(test_corr_data['Age'], test_corr_data['Error'], color='blue')

        plt.title('True age vs Error')
        plt.xlabel('True age')
        plt.ylabel('Error')
        plt.grid(True)
        plt.show()
        
        # Plotting true vs predicted values with smoothed curve and ElasticNet predictions
        plt.figure(figsize=(10, 6))
        plt.scatter(test_corr_data['Age'], test_corr_data['Predicted Age'], color='blue', label='True vs Predicted')
        plt.plot(test_corr_data_sorted['Age'], smoothed_values, color='red', label='Smoothed Curve')

        plt.title('True vs Predicted Values with Smoothed Curve')
        plt.xlabel('True Age')
        plt.ylabel('Predicted Age')
        plt.legend()
        plt.grid(True)
        plt.show()

    return test_corr_data


def correlation_analysis(test_corr_data, metadata, override=False):
    colnames_new = []
    p_vec = []
    pval_vec = []
    
    if override:
        errors = test_corr_data
    else:
        errors = test_corr_data['Individual age gaps']

    for i in np.arange(0, len(metadata.columns)):
        predictor = metadata.iloc[:, i]
        non_nan_idx = ~predictor.isna()
        predictor = predictor[non_nan_idx]
        errors_filtered = errors[non_nan_idx]

        if len(errors_filtered) >= 2 and len(predictor) >= 2:
            uniques = predictor.unique()

            # check if all values are numeric, then calculate pearson correlation
            if all((np.issubdtype(type(x), np.integer) or np.issubdtype(type(x), np.floating)) for x in uniques):

                pearson_corr, pval = pearsonr(errors_filtered, predictor)
                if not np.isnan(pearson_corr):
                    p_vec.append(pearson_corr)
                    pval_vec.append(pval)
                    colnames_new.append(metadata.columns[i])


    corr_result = pd.DataFrame({'Predictor':colnames_new, 'Pearson Correlation':p_vec, 'P-value': pval_vec})
    corr_result = corr_result.iloc[np.abs(corr_result['Pearson Correlation']).argsort()[::-1]]

    # FDR correction
    reject, corrected_pvals, _, _ = multitest.multipletests(corr_result['P-value'], method='fdr_bh')
    corr_result['Corrected P-value'] = corrected_pvals

    # Test null hypothesis
    corr_result['Reject Null Hypothesis'] = reject

    # Sort by corrected p-value
    corr_result = corr_result.iloc[corr_result['Corrected P-value'].argsort()]
    corr_result.reset_index(drop=True, inplace=True)

    return corr_result

def get_metadata_average(R1_Meta, R2_Meta, rows):
    R1_Metadata = R1_Meta.copy()
    R2_Metadata = R2_Meta.copy()

    R1_Metadata.columns = [col + '_R1'if col != 'visit_id' else col for col in R1_Metadata.columns]
    R2_Metadata.columns = [col + '_R2'if col != 'visit_id' else col for col in R2_Metadata.columns]

    all_metadata = pd.merge(R1_Metadata, R2_Metadata, left_index=True, right_index=True)

    # tmp_metadata = all_metadata.iloc[y_test_avg.index].reset_index()
    tmp_metadata = all_metadata.loc[rows]
    metadata_average = pd.DataFrame()
    metadata_average['Carotid IMT average'] = (tmp_metadata['Carotid IMT R2_R2'] + tmp_metadata['Carotid IMT R1_R2'])/2
    metadata_average['Femoral IMT average'] = (tmp_metadata['Femoral IMT R2_R2'] + tmp_metadata['Femoral IMT R1_R2'])/2
    metadata_average['BMI averaged (kg/m2)'] = tmp_metadata['BMI averaged (kg/m2)_R2']
    metadata_average['Systolic BP average (mmHg)'] = tmp_metadata['Systolic BP averaged (mmHg)_R2']
    metadata_average['Diastolic BP average (mmHg)'] = tmp_metadata['Diastolic BP averaged (mmHg)_R2']
    metadata_average['Pulse Pressure average (mmHg)'] = tmp_metadata['Pulse pressure averaged (mmHg)_R2']
    metadata_average['Total cholesterol average'] = tmp_metadata['Total cholesterol averaged_R2']
    metadata_average['HDL cholesterol average'] = tmp_metadata['HDL cholesterol averaged_R2']
    metadata_average['LDL cholesterol average'] = tmp_metadata['LDL cholesterol averaged_R2']
    metadata_average['Triglycerides average (mg/dl)'] = (tmp_metadata['Triglycerides (mg/dl; R2)_R2'] + tmp_metadata['Triglycerides (mg/dl; R1)_R2'])/2
    metadata_average['High sensitive CRP average (mg/dL)'] = tmp_metadata['High sensitive CRP averaged (mg/dL)_R2']
    metadata_average['IL6 average'] = (tmp_metadata['IL6 R2_R2'] + tmp_metadata['IL6 R1_R2'])/2
    metadata_average['White blood cells average'] = (tmp_metadata['WBC R2_R2'] + tmp_metadata['Whitebloodcellcount10e3microl_R1'])/2
    metadata_average['Fibrinogen average'] = (tmp_metadata['Fibrinogen R2_R2'] + tmp_metadata['Fibrinogen R1_R2'])/2
    metadata_average['Homocystein average'] = (tmp_metadata['Homocystein R2_R2'] + tmp_metadata['Homocystein R1_R2'])/2
    metadata_average['Glycemia average'] = tmp_metadata['Glycemia averaged_R2']
    metadata_average['Creatinin average (mg/dl)'] = (tmp_metadata['creat R2_R2'] + tmp_metadata['Creatinin_mgdl_R1'])/2
    metadata_average['Uric acid average (mg/dL)'] = tmp_metadata['Uric acid averaged (mg/dL)_R2']
    metadata_average['Pack years nicotine average (year)'] = (tmp_metadata['Pack years nicotine R1 (year)_R2'] + tmp_metadata['Pack years nicotine R2 (year)_R2'])/2

    return metadata_average



def get_metadata_delta(R1_Meta, R2_Meta, rows):
    R1_Metadata = R1_Meta.copy()
    R2_Metadata = R2_Meta.copy()

    R1_Metadata.columns = [col + '_R1'if col != 'visit_id' else col for col in R1_Metadata.columns]
    R2_Metadata.columns = [col + '_R2'if col != 'visit_id' else col for col in R2_Metadata.columns]

    all_metadata = pd.merge(R1_Metadata, R2_Metadata, left_index=True, right_index=True)

    tmp_metadata = all_metadata.loc[rows]
    metadata_delta = pd.DataFrame()
    metadata_delta['Carotid IMT delta'] = tmp_metadata['Carotid IMT delta_R2']
    metadata_delta['Femoral IMT delta'] = tmp_metadata['Femoral IMT delta_R2']
    metadata_delta['BMI delta (kg/m²)'] = tmp_metadata['BMI delta (kg/m2)_R2']
    metadata_delta['Systolic BP delta (mmHg)'] = tmp_metadata['Systolic BP delta (mmHg)_R2']
    metadata_delta['Diastolic BP delta (mmHg)'] = tmp_metadata['Diastolic BP delta (mmHg)_R2']
    metadata_delta['Pulse Pressure delta (mmHg)'] = tmp_metadata['Pulse Pressure delta (mmHg)_R2']
    metadata_delta['Total cholesterol delta'] = tmp_metadata['Total cholesterol delta_R2']
    metadata_delta['HDL cholesterol delta'] = tmp_metadata['HDL cholesterol delta_R2']
    metadata_delta['LDL cholesterol delta'] = tmp_metadata['LDL cholesterol delta_R2']
    metadata_delta['Triglycerides delta (mg/dl)'] = tmp_metadata['Triglycerides (mg/dl; R2)_R2'] - tmp_metadata['Triglycerides (mg/dl; R1)_R2']
    metadata_delta['High sensitive CRP delta (mg/dL)'] = tmp_metadata['High sensitive CRP delta (mg/dL)_R2']
    metadata_delta['IL6 delta'] = tmp_metadata['IL6 R2_R2'] - tmp_metadata['IL6 R1_R2']
    metadata_delta['White blood cells delta'] = tmp_metadata['WBC R2_R2'] - tmp_metadata['Whitebloodcellcount10e3microl_R1']
    metadata_delta['Fibrinogen delta'] = tmp_metadata['Fibrinogen R2_R2'] - tmp_metadata['Fibrinogen R1_R2']
    metadata_delta['Homocystein delta'] = tmp_metadata['Homocystein R2_R2'] - tmp_metadata['Homocystein R1_R2']
    metadata_delta['Glycemia delta'] = tmp_metadata['Glycemia delta_R2']
    metadata_delta['Creatinin delta (mg/dl)'] = tmp_metadata['creat R2_R2'] - tmp_metadata['Creatinin_mgdl_R1']
    metadata_delta['Uric acid delta (mg/dL)'] = tmp_metadata['Uric acid delta (mg/dL)_R2']
    metadata_delta['Pack years nicotine delta (year)'] = tmp_metadata['Pack years nicotine delta (year)_R2']
    
    return metadata_delta

def get_mixed_metadata(R1_Meta, R2_Meta, matched_data, chosen_rounds):
    R1_Metadata = R1_Meta.copy()
    R2_Metadata = R2_Meta.copy()

    R1_Metadata.columns = [col + '_R1'if col != 'visit_id' else col for col in R1_Metadata.columns]
    R2_Metadata.columns = [col + '_R2'if col != 'visit_id' else col for col in R2_Metadata.columns]

    all_metadata = pd.merge(R1_Metadata, R2_Metadata, left_index=True, right_index=True)

    metadata_individuals = pd.DataFrame(columns=['Carotid IMT', 'Femoral IMT', 'BMI (kg/m²)', 'Systolic BP (mmHg)', 'Diastolic BP (mmHg)',
                                             'Pulse Pressure (mmHg)', 'Total cholesterol', 'HDL cholesterol', 'LDL cholesterol',
                                             'Triglycerides (mg/dl)', 'High sensitive CRP (mg/dL)', 'IL6', 'White blood cells', 'Fibrinogen',
                                             'Homocystein', 'Glycemia', 'Creatinin (mg/dl)', 'Uric acid (mg/dL)', 'Pack years nicotine (year)', 'Smoker'])

    columns_r1 = ['Carotid IMT R1_R2', 'Femoral IMT R1_R2', 'BMI R1 (kg/m2)_R2', 'Systolic BP R1 (mmhg)_R2', 'Diastolic BP R1 (mmHg)_R2',
                                                'Pulse Pressure R1 (mmHg)_R2', 'Total cholesterol R1_R2', 'HDL cholesterol R1_R2', 'LDL cholesterol R1_R2',
                                                'Triglycerides (mg/dl; R1)_R2', 'High sensitive CRP R1 (mg/dL)_R2', 'IL6 R1_R2', 'Whitebloodcellcount10e3microl_R1', 
                                                'Fibrinogen R1_R2', 'Homocystein R1_R2', 'Glycemia R1_R2', 'Creatinin_mgdl_R1', 'Uric acid R1 (mg/dL)_R2',
                                                'Pack years nicotine R1 (year)_R2', 'Current Smoker R1_R2']

    columns_r2 = ['Carotid IMT R2_R2', 'Femoral IMT R2_R2', 'BMI R2 (kg/m2)_R2', 'Systolic BP R2 (mmHg)_R2', 'Diastolic BP R2 (mmHg)_R2',
                                                'Pulse Pressure R2 (mmHg)_R2', 'Total cholesterol R2_R2', 'HDL cholesterol R2_R2', 'LDL cholesterol R2_R2',
                                                'Triglycerides (mg/dl; R2)_R2', 'High sensitive CRP R2 (mg/dL)_R2', 'IL6 R2_R2', 'WBC R2_R2', 
                                                'Fibrinogen R2_R2', 'Homocystein R2_R2', 'Glycemia R2_R2', 'creat R2_R2', 'Uric acid R2 (mg/dL)_R2',
                                                'Pack years nicotine R2 (year)_R2', 'Current Smoker R2_R2']

    for i, visit_id in enumerate(matched_data.index):
        chosen_round = chosen_rounds[i]
        if chosen_round == 1:
            metadata_individuals.loc[visit_id] = all_metadata.loc[visit_id][columns_r1].values
        else:
            metadata_individuals.loc[visit_id] = all_metadata.loc[visit_id][columns_r2].values

    return metadata_individuals