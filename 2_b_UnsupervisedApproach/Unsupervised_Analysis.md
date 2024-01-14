# Unsupervised Analysis using Dynamic Tree Clustering and evaluating cluster quality through Silhouette Score


This repository contains a Jupyter Notebook for performing unsupervised analysis on a dataset, primarily focusing on clustering quality evaluation using the Silhouette Score. The notebook includes a comprehensive set of Python scripts for data processing, clustering, and result visualization.

## Getting Started

To run this Jupyter Notebook, create a new Conda environment using the `requirements.txt` file located in the same GitHub folder. This file lists all the necessary Python packages and their versions.

### Prerequisites

- Anaconda or Miniconda
- Python 3.x

### Installation

1. **Download the repository:**

   Download the ZIP file of the repository and extract it, or copy the files to your working directory.

2. **Create a new Conda environment:**

   Navigate to the directory containing the `requirements.txt` file and create a new Conda environment:

   ```bash
   conda create --name <env-name> --file requirements.txt

### Activate the Conda environment:

Activate the new environment:

```bash
conda activate <env-name>


## Running the Notebook

After installing the prerequisites, you can run the Jupyter Notebook.

1. **Launch Jupyter Notebook:**

    ```bash
    jupyter notebook
    ```

2. **Open the `Unsupervised_Analysis.ipynb` file in the Jupyter Notebook interface.**

3. **Run the cells sequentially to perform the analysis.**

## File Structure

- `Unsupervised_Analysis.ipynb` - Main Jupyter Notebook containing the analysis scripts.
- `requirements.txt` - Lists all the Python packages required to run the notebook.
- `/data` - Directory containing the datasets and other data files used in the notebook.

## Notebook Overview

The Jupyter Notebook is divided into various sections, each performing specific tasks:

1. **Import Modules**: Importing necessary Python libraries and modules.
2. **Setting Base Directory**: Setting up the working directory for file operations.
3. **Reading Files**: Loading the dataset and other necessary files.
4. **Data Preprocessing**: Preprocessing the data for analysis.
5. **Clustering**: Performing clustering on the preprocessed data.
6. **Silhouette Score Evaluation**: Calculating and analyzing the Silhouette Score to evaluate the quality of clustering.
7. **Result Visualization**: Visualizing the results of the analysis through plots and graphs.

## Contributing

Contributions to this project are welcome. Please ensure to update tests as appropriate.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
