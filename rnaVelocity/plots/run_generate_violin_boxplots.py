import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu

#VLNPLT

# Define a function to extract relevant data from each AnnData object
def extract_data_from_h5ad(file_path, index_str = "1"):
    adata = sc.read_h5ad(file_path)
    
    data = {
    'index2': adata.obs.index,
    'Regions': adata.obs['Regions'].values,
    'Condition': adata.obs['Condition'].values,
    'ConditionW': [cond + index_str for cond in adata.obs['Condition']],
    'pxl_row_in_fullres': adata.obs['pxl_row_in_fullres'].values,
    'pxl_col_in_fullres': adata.obs['pxl_col_in_fullres'].values,
    'sample_batch': adata.obs['sample_batch'].values,
    'n_counts': adata.obs['n_counts'].values,
    'Sample': adata.obs['Sample'].values,
    'U_US_ratio': adata.obs['U_US_ratio'].values,
    'latent_time': adata.obs['latent_time'].values,
    'velocity_length': adata.obs['velocity_length'].values,
    'velocity_pseudotime': adata.obs['velocity_pseudotime'].values,
    }
    
    df = pd.DataFrame(data)
    #df.set_index('index', inplace=True)
    #df = adata.obs[['Regions', 'Condition', 'Sample', 'U_US_ratio', 'latent_time', 'velocity_pseudotime', 'velocity_length']]
    print("Processed: " + file_path )
    return df


def save_dataframe_to_disk(df, file_path):
    """
    Save the concatenated DataFrame to disk.

    Args:
    df (pd.DataFrame): DataFrame to save.
    file_path (str): Path to save the DataFrame.
    """
    if os.path.exists(file_path):
        return
    df.to_csv(file_path, index=False)
    print(f"DataFrame saved to {file_path}")


def read_dataframe_from_disk(file_path):
    """
    Read the concatenated DataFrame from disk.

    Args:
    file_path (str): Path to the saved DataFrame.

    Returns:
    pd.DataFrame: DataFrame read from disk.
    """
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        print(f"DataFrame read from {file_path}")
        return df
    raise Exception(f"Failed to load {file_path}")


def wilcoxon_rank_sum_test(df, unit, dataset_type="F", within = False):
    """
    Perform Wilcoxon Rank-Sum Test (Mann-Whitney U Test) for each region.

    Args:
    df (pd.DataFrame): DataFrame containing the data to test.

    Returns:
    pd.DataFrame: DataFrame containing the test results.
    """

    condition_values = {
            'F': {'0':"CF", '1':'EcigF'},
            'M':{'0':"CM", '1':'EcigM'},
            }

    condition_valuesW = {
            'F': {'0':"CF1", '1':'CF2'},
            'M':{'0':"CM1", '1':'CM2'},
            }


    results = []
    listc1 = [region + "_CM1" for region in df['Regions'].unique()]
    listc2 = [region + "_CM2" for region in df['Regions'].unique()]

    plot_df = pd.DataFrame()

    if not within:
        for region in df['Regions'].unique():
            control_data = df[(df['Regions'] == region) & (df['Condition'] == condition_values[dataset_type]['0'])][unit].values
            treated_data = df[(df['Regions'] == region) & (df['Condition'] == condition_values[dataset_type]['1'])][unit].values
            if len(control_data) > 0 and len(treated_data) > 0:
                stat, p_value = mannwhitneyu(control_data, treated_data, alternative='two-sided')
                results.append({'Region': region, 'U-statistic': stat, 'p-value': p_value})
    else:
        for region in df['Regions'].unique():
            control_data = df[(df['Regions'] == region) & (df['ConditionW'] == condition_valuesW[dataset_type]['0'])][unit].values
            treated_data = df[(df['Regions'] == region) & (df['ConditionW'] == condition_valuesW[dataset_type]['1'])][unit].values

            if len(control_data) > 0 and len(treated_data) > 0:
                print(region)
                print(control_data)
                print(region)
                print(treated_data)

                stat, p_value = mannwhitneyu(control_data, treated_data, alternative='two-sided')
                results.append({'Region': region, 'U-statistic': stat, 'p-value': p_value})
        
    #print(results)
    return pd.DataFrame(results)


def print_wilcoxon_results(results_df, unit, dataset="FA"):
    """
    Print the Wilcoxon Rank-Sum test results for each region.

    Args:
    results_df (pd.DataFrame): DataFrame containing the test results.
    """
    sorted_results_df = results_df.sort_values(by='p-value')
    for _, row in sorted_results_df.iterrows():
        print(f"{dataset},{unit},{row['Region']},{row['U-statistic']:.3f},{row['p-value']:.3e},")


def plot_distributions_for_regions(df, unit="latent_time", dataset="MA", output_dir='distributions'):
    """
    Plot distributions for each Region, comparing distributions across Conditions using violin plots.

    Args:
    df (pd.DataFrame): DataFrame containing the data to plot.
    unit (str): The variable to plot on the y-axis.
    dataset (str): The name of the dataset for the title and filename.
    output_dir (str): Directory where the plots will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.figure(figsize=(18, 12))

    # Violin plot with Condition as hue
    sns.violinplot(data=df, x='Regions', y=unit, hue='Condition', split=True, inner='quart', palette='muted', legend=False)

    # Strip plot with Sample as hue, using a different palette for Samples
    num_samples = len(df['Sample'].unique())
    palette = sns.color_palette("hsv", num_samples)
    strip_plot = sns.stripplot(data=df, x='Regions', y=unit, hue='Sample', dodge=True, jitter=True, palette=palette, alpha=0.7, size=3, marker='o')

    # Add text labels next to each distribution
    for region in df['Regions'].unique():
        region_data = df[df['Regions'] == region]
        for sample in region_data['Sample'].unique():
            sample_data = region_data[region_data['Sample'] == sample]
            median_value = sample_data[unit].median()
            x_position = df['Regions'].unique().tolist().index(region)
            strip_plot.text(x=x_position, y=median_value, s=sample, color='black', ha='center')

    # Customize the plot
    plt.title(f"Distribution of {unit} in {dataset} for each Region and Sample")
    plt.xlabel('Region')
    plt.ylabel(unit)
    plt.xticks(rotation=45)

    # Adjust the legend to show both Condition and Sample
    handles, labels = plt.gca().get_legend_handles_labels()
    n_conditions = len(df['Condition'].unique())
    n_samples = len(df['Sample'].unique())

    # Remove duplicate legend entries for Samples
    sample_handles = handles[-n_samples:]
    sample_labels = labels[-n_samples:]
    plt.legend(sample_handles, sample_labels, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    #plt.legend(handles, labels, title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Save the plot
    output_path = os.path.join(output_dir, f"VLNPLT_Jan3_{dataset}_{unit}.pdf")
    plt.savefig(output_path, format='pdf', dpi=1200)
    plt.close()
    plot_boxplots_for_regions(df, unit, dataset, output_dir)


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_distributions_for_regions2(
    df_original, 
    unit="latent_time", 
    dataset="MA", 
    output_dir='distributions', 
    use_stripplot=False
):
    """
    Plot distributions for each Region, comparing distributions across Conditions using violin plots.
    Optionally add a strip plot on top by setting use_stripplot=True.

    Args:
        df (pd.DataFrame): DataFrame containing the data to plot. Must have columns 'Regions', 'Condition', 'Sample'.
        unit (str): The variable to plot on the y-axis.
        dataset (str): The name of the dataset for the title and filename.
        output_dir (str): Directory where the plots will be saved.
        use_stripplot (bool): Whether to add a stripplot (showing individual points) on top of the violin plot.
    """

    unique_regions = []
    if "A" in dataset:
        unique_regions = ['LV','ACC','Septal','Pir','AIC','NAc','CPu','cc','Layer5/6','Layer4','Layer2/3','Layer1']
    else:
        unique_regions = ['LV','Amy','HY','TH','HP','RSP','CPu','cc','Layer5/6','Layer4','Layer2/3','Layer1']

    df = df_original[df_original['Regions'].isin(unique_regions)]

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.figure(figsize=(18, 12))

    # If we are using a strip plot, we don't want the violinplot to manage a legend for 'Condition',
    # because we'll manage legends manually later. If not using stripplot, we can let the violin
    # plot create a legend.
    violin_legend_setting = not use_stripplot

    # Violin plot with Condition as hue
    sns.violinplot(
        data=df, 
        x='Regions', 
        y=unit, 
        hue='Condition', 
        split=True, 
        inner='quart', 
        palette='muted', 
        legend=violin_legend_setting,
        order=unique_regions,
    )

    # Optionally add a strip plot for Samples
    if use_stripplot:
        num_samples = df['Sample'].nunique()
        palette = sns.color_palette("hsv", num_samples)

        strip_plot = sns.stripplot(
            data=df,
            x='Regions',
            y=unit,
            hue='Sample',
            dodge=True,
            jitter=True,
            palette=palette,
            alpha=0.7,
            size=3,
            marker='o'
        )

        # Add text labels next to each distribution (showing the Sample name at the median)
        for region in unique_regions:
            region_data = df[df['Regions'] == region]
            x_position = unique_regions.index(region)
            for sample in region_data['Sample'].unique():
                sample_data = region_data[region_data['Sample'] == sample]
                median_value = sample_data[unit].median()
                strip_plot.text(
                    x=x_position,
                    y=median_value,
                    s=sample,
                    color='black',
                    ha='center'
                )

    # Customize the plot
    plt.title(f"Distribution of {unit} in {dataset} for each Region")
    plt.xlabel('Region')
    plt.ylabel(unit)
    plt.xticks(rotation=45)

    # Legend handling
    if use_stripplot:
        # We have two different hue encodings: Condition (for violin) and Sample (for stripplot).
        # Extract the legend handles and labels from the current axes
        handles, labels = plt.gca().get_legend_handles_labels()
        n_samples = df['Sample'].nunique()

        # The last n_samples handles/labels belong to 'Sample' 
        # (assuming Seaborn layered them in that order).
        sample_handles = handles[-n_samples:]
        sample_labels = labels[-n_samples:]

        # Everything else before those belongs to 'Condition'.
        # But we may or may not show them here. For simplicity, we can
        # just show the sample legend. If you also want the condition legend,
        # you can handle them accordingly.
        plt.legend(
            sample_handles, 
            sample_labels, 
            title='Sample',
            bbox_to_anchor=(1.05, 1), 
            loc='upper left'
        )

    else:
        # If not using stripplot, the violinplotâ€™s legend shows Condition (unless the data is two categories only).
        # Seaborn should have automatically created it if legend=violin_legend_setting is True.
        plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Save the plot
    output_path = os.path.join(output_dir, f"VLNPLT_Jan7_{dataset}_{unit}.pdf")
    plt.savefig(output_path, format='pdf', dpi=1200)
    plt.close()

    # Call boxplot function (make sure this is defined or imported in your environment)
    plot_boxplots_for_regions(df, unit, dataset, output_dir, unique_regions)


def plot_distributions_for_regions_legacy(df, unit="latent_time", dataset="MA", output_dir='distributions'):
    """
    Plot distributions for each Region, comparing distributions across Conditions using violin plots.

    Args:
    df (pd.DataFrame): DataFrame containing the data to plot.
    unit (str): The variable to plot on the y-axis.
    dataset (str): The name of the dataset for the title and filename.
    output_dir (str): Directory where the plots will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    plt.figure(figsize=(18, 12))
    sns.violinplot(data=df, x='Regions', y=unit, hue='Condition', split=True, inner='quart', palette='muted')

    # Add strip plot for each Sample, differentiated by Condition
    sns.stripplot(data=df, x='Regions', y=unit, hue='Sample', dodge=True, jitter=True, palette='dark:.3', alpha=0.5, size=3, marker='o')

    # Customize the plot
    plt.title(f"Distribution of {unit} in {dataset} for each Region and Sample")
    plt.xlabel('Region')
    plt.ylabel(unit)
    plt.xticks(rotation=45)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    n = len(df['Condition'].unique())
    plt.legend(handles[:n], labels[:n], title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the plot
    output_path = os.path.join(output_dir, f"VLNPLT_Jan3_{dataset}_{unit}.pdf")
    plt.savefig(output_path, format='pdf', dpi=1200)


# Define a function to plot boxplots for each Region
def plot_boxplots_for_regions(df, unit = "latent_time", dataset="MA", output_dir='boxplots', unique_regions = []):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.figure(figsize=(14, 10))
    sns.boxplot(data=df, x='Regions', y=unit, hue='Condition', order=unique_regions)
    # Customize the plot
    plt.title("Boxplot of " + unit + " in " + dataset + " for each Region")
    plt.xlabel('Region')
    plt.ylabel(unit)
    plt.xticks(rotation=45)
    plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"BXPLT_Jan7_{dataset}_{unit}.pdf")
    plt.savefig(output_path, format='pdf', dpi=1200)
    # Show the plot
    plt.show()


# Example usage
file_ma_paths = [
        './VersionJuly1/velo_ma_anndata_saved/P94A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P95A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P100A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P80A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P83A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P112A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P81A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P113A1.h5ad',
        './VersionJuly1/velo_ma_anndata_saved/P81A2.h5ad',
        ]  # Replace with your file paths

file_deconv_ma_paths = [
        './e_i_ratio_d/P94A1.decov.csv',
        './e_i_ratio_d/P95A1.decov.csv',
        './e_i_ratio_d/P100A1.decov.csv',
        './e_i_ratio_d/P80A1.decov.csv',
        './e_i_ratio_d/P83A1.decov.csv',
        './e_i_ratio_d/P112A1.decov.csv',
        './e_i_ratio_d/P81A1.decov.csv',
        './e_i_ratio_d/P113A1.decov.csv',
        './e_i_ratio_d/P81A2.decov.csv',
        ]  # Replace with your file paths


file_mp_paths = [
    './VersionJuly1/velo_mp_anndata_saved/P94P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P95P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P100P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P80P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P83P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P112P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P81P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P113P1.h5ad',
    './VersionJuly1/velo_mp_anndata_saved/P81P2.h5ad',
    ]  # Replace with your file paths # Initialize an empty DataFrame to collect all data

# Example usage
file_fa_paths = [
        './VersionJuly1/velo_fa_anndata_saved/P96A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P93A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P101A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P108A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P111A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P74A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P77A1.h5ad',
        './VersionJuly1/velo_fa_anndata_saved/P82A1.h5ad',
        #'./velo_fa_anndata_saved/P95A1.h5ad'
        ]  # Replace with your file paths


file_deconv_fa_paths = [
        './e_i_ratio_d/P96A1.decov.csv',
        './e_i_ratio_d/P93A1.decov.csv',
        './e_i_ratio_d/P101A1.decov.csv',
        './e_i_ratio_d/P108A1.decov.csv',
        './e_i_ratio_d/P111A1.decov.csv',
        './e_i_ratio_d/P74A1.decov.csv',
        './e_i_ratio_d/P77A1.decov.csv',
        './e_i_ratio_d/P82A1.decov.csv'
        ]

file_fp_paths = [
        './VersionJuly1/velo_fp_anndata_saved/P101P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P108P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P111P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P74P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P77P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P82P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P93P1.h5ad',
        './VersionJuly1/velo_fp_anndata_saved/P96P1.h5ad'
    ]

import sys
import numpy as np


def create_pd_from_files_extended(file_paths, file_deconv_paths, file_path = "collected_saved_anndata.dfp"):
    tmp = pd.DataFrame()
    # Extract data from each file and append it to the DataFrame

    try:
        tmp = read_dataframe_from_disk(file_path)
        return tmp
    except:
        tmp = pd.DataFrame()
        tmp2 = pd.DataFrame()
        for i in range(len(file_paths)):
            df1 = extract_data_from_h5ad(file_paths[i], str(i))
            df2 = pd.read_csv(file_deconv_paths[i])
            df2['e_i_ratio'] = df2["excit_N"] / df2["inhib_N"]
            df2['log_e_i_ratio'] = -1*np.log(df2['e_i_ratio']) 
            df2['index2'] = df2['ind'].str[4:16]
            df2['sample_batch'] =  df2['ind'].str[-6:]
            df = pd.merge(df1, df2, on=['index2', 'sample_batch'])
            tmp = pd.concat([tmp, df], ignore_index=True)

        save_dataframe_to_disk(tmp, file_path)
        return tmp


def create_pd_from_files(file_paths, file_path = "collected_saved_anndata.dfp"):
    tmp = pd.DataFrame()
    # Extract data from each file and append it to the DataFrame
    
    try: 
        tmp = read_dataframe_from_disk(file_path)
        return tmp
    except:
        tmp = pd.DataFrame()
        for (idx, filename) in enumerate(file_paths):
            df = extract_data_from_h5ad(filename, str(idx))
            tmp = pd.concat([tmp, df], ignore_index=True)
        
        save_dataframe_to_disk(tmp, file_path)
        return tmp

import sys


#units = ["latent_time","velocity_pseudotime","velocity_length","U_US_ratio", "log_e_i_ratio"]
units = ["latent_time","velocity_pseudotime","velocity_length","U_US_ratio"]

if True:
    all_data_ma_df = create_pd_from_files_extended(file_ma_paths, file_deconv_ma_paths, "ma_saved_merged_ext_annData_w.dfp")
# Plot the boxplots
    for unit in units:
        plot_distributions_for_regions2(all_data_ma_df, unit, "MA")

# Perform Wilcoxon Rank-Sum test (Latent-time)
    for unit in units:
        wilcoxon_results_unit_df = wilcoxon_rank_sum_test(all_data_ma_df, unit, "M")
        #print(f"Wilcoxon Rank-Sum test results (Latent-time):")
        print_wilcoxon_results(wilcoxon_results_unit_df,unit, "MA")


if True:
    all_data_mp_df = create_pd_from_files(file_mp_paths, "mp_saved_merged_annData_w.dfp")
    for unit in units:
        plot_distributions_for_regions2(all_data_mp_df, unit, "MP")

# Perform Wilcoxon Rank-Sum test (Latent-time)
    for unit in units:
        wilcoxon_results_unit_df = wilcoxon_rank_sum_test(all_data_mp_df, unit, "M")
        #print(f"Wilcoxon Rank-Sum test results (Latent-time):")
        print_wilcoxon_results(wilcoxon_results_unit_df, unit, "MP")


if True:
    all_data_fa_df = create_pd_from_files_extended(file_fa_paths, file_deconv_fa_paths, "fa_saved_merged_ext_annData_w.dfp")
# Plot the boxplots
    for unit in units:
        plot_distributions_for_regions2(all_data_fa_df, unit, "FA")

# Perform Wilcoxon Rank-Sum test (Latent-time)
    for unit in units:
        wilcoxon_results_unit_df = wilcoxon_rank_sum_test(all_data_fa_df, unit, "F")
        #print(f"Wilcoxon Rank-Sum test results (Latent-time):")
        print_wilcoxon_results(wilcoxon_results_unit_df, unit, "FA")


if True:
    all_data_fp_df = create_pd_from_files(file_fp_paths, "fp_saved_merged_annData_w.dfp")
# Plot the boxplots
    for unit in units:
        plot_distributions_for_regions2(all_data_fp_df, unit, "FP")

# Perform Wilcoxon Rank-Sum test (Latent-time)
    for unit in units:
        wilcoxon_results_unit_df = wilcoxon_rank_sum_test(all_data_fp_df, unit, "F")
        #print(f"Wilcoxon Rank-Sum test results (Latent-time):")
        print_wilcoxon_results(wilcoxon_results_unit_df, unit, "FP")

