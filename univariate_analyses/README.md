
A repository for the correlation/machine learning analyses on the Dynamic
Brain and Mind Lab's EEG-fMRI data.


# Installation

This project requires both **R** and **Python**, with Python libraries
managed through **Miniconda**.

## Step 1: Install R

### Option A: Homebrew (recommended)

If you use [Homebrew](https://brew.sh/), install R with

`brew install --cask r`

This installs the latest stable version of R system wide. 

### Option B: via CRAN

If you don't use Homebrew:

1. Visit https://cran.r-project.org
2. Click "Download R for macOS"
3. Download and install the `.pkg` for the latest release

Once installed, verify:

```bash
R --version
```
You should see something like `R version 4.x.x (YYYY-MM-DD)`.

## Step 2 (Optional): Install RStudio

Rstudio provies a user-friendly IDE for R, similar to the built in IDE for
Matlab.

Download from [here](https://posit.co/download/rstudio).

## Step 3: Install Miniconda

Miniconda will manage your Python environment (and will also install
Python). 

1. Download the latest **Miniconda for macOS** from
   [here](https://docs.conda.io/en/latest/miniconda.html). This will
   download a `.sh` script that you will need to navigate to in a terminal
   and run with something like `bash Miniconda3-latest-MacOSX-arm64.sh`.

   You could alternatively download *and* install via terminal with:

   ```bash
   curl -fsSLo ~/Downloads/Miniconda3-latest-MacOSX-arm64.sh \
     https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
   bash ~/Downloads/Miniconda3-latest-MacOSX-arm64.sh
   ```

2. Follow the prompts and allow it to initialize your shell (adds lines to
   `.zshrc` or `.bashrc`).

3. Close and reopen your terminal, then verify:

```bash
conda --version
```
You should see something like `conda 24.x.x`.


## Step 4: Create the project environment

### Configure Python environment

Once Miniconda is installed and working, you can create the project's
Python environment by configuring and activating the Python environment for
this project (which will install all needed code libraries).

```bash
cd ~ # Or wherever you want to clone the repo
git clone https://github.com/dbraun31/eeg-fmri-ml.git
cd eeg-fmri-ml
conda env create -f environment.yml
conda activate eeg-fmri
```
You can test whether everything installed correctly in two ways:

```bash
python --version
```
You should see `python 3.11.x`.

```bash
python test_python_imports.py
```
This script will output a message indicating whether all libraries can be
imported successfully.

*Note* that whenever you open a new terminal (and see `(base)` in your
command line), you will need to run `conda activate eeg-fmri` to 'turn on'
this environment. You will get package import errors if you forget this
step.


### Configure R environment

*Note:* I did not choose to configure a dedicated R virtual environment for this
project. Here, I'm showing you how to install all packages in your global
environment (I have far fewer package conflicts in R vs. Python; but you
could use an R env manager if you prefer).

To install R packages, run:

```bash
Rscript environment.r
```

Unless you are very lucky, this script is unlikely to install all packages
successfully on your first try. The following two packages in particular
are very finicky and likely to require dependency hunting and
troubleshooting: `arrow`, `eegUtils`. The second package (`eegUtils`) can
be extremely painful, and it's only necessary for making EEG topo plots in R.
You can skip installing `eegUtils` if you don't care about topo plots, or
follow the troubleshooting instructions in `eeg_utils_troubleshooting.md`
to go down the rabbit hole (we can work on this one together).

To verify whether packages have installed successfully, run:

```bash
Rscript test_r_imports.r
```

