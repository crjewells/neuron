# neuron
Catalogue of research of dendritic input-output nonlinearity in pyramidal neurons using NEURON module


# Environment Setup

Before starting to code, there are a few steps we need to take to setup our python environment in order for NUERON to run properly.

## Install Brew

If you haven't already, install brew by running the following in your terminal. This is helpful installer.

`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`

## Virtual Environment

We will want to download a virtual environment software in order to bundle everything we need to run in a safe, isolated environment. Make sure you have Python downloaded, which offers virtual environment capabilites.

Once the application has downloaded, we can set up a virtual environment in command line. This creates a new space where we can download all the packages needed for our programming.

In terminal:

Navigate to whichever directory you want to setup your environment in. This doesn't really matter, but it's probably better to put in the same directory as wherever you plan to store your python files.

`python3 -m venv <env-name> # creates a new environment`

Then activate your environment:

`source <env-name>/bin/activate`

You can deactivate your environment anytime by using:

`deactivate`

## Installing NEURON

You'll want to refer to this [guide](https://nrn.readthedocs.io/en/8.2.6/install/install_instructions.html#mac-os) if you get stuck, but installation is pretty straight forward.

We first will install the NEURON package in our virtual environment, so make sure it is activated.

`pip3 install neuron # installs neuron package`

If that doesn't work, you can also try pip install. The system will tell you if you need to update your pip/pip3 installer. It may also be likely that you have to run this outside of your environment (install system wide), but I'm not to sure since I haven't done it in a while.

## Installing Command Line Tools and XQuartz

NEURON needs more packages to run, so you'll need to download these as well.

Install command line tools:

`xcode-select --install`

You can download XQuartz [here](https://www.xquartz.org)`

## Installing Jupyter Notebook

You can run your code using Jupyeter notebook which is helpful when running and testing chunks of code

Install jupyter using (in your virtual environment):

`pip install jupyter`

Install python3 kernel using:

`pip install ipykernel`

# Importing Mods

You will often want to download mods, like NMDA synapes or other channels, that are not native to neuron. 

1. Place all .mod files DIRECTLY into the folder that contains the Python files you would like to run

2. From terminal, activate your virtual environment

3. Navigate to the folder that contains your mods through terminal. You can do so by the `cd` command
    1. Using `cd /XX/XX` will direct you to an absolute path
    2. Using `ls`, you can find the list of directories within the folder you are in, and then use `cd <directory-name>` to access the directory

4. Once you are in the folder, type `sudo nrnivmodl` to compile the mods

5. From here, you should be able to run Jupyter directly; as long as the Python files are in the same folder that has “arm64”, you should be good!


# Using CARC

## Setting Up CARC Environment

When you start complex simulations, your run time will start to get ridiculous. Towards the end, I began running my simulations through USC's advanced computers, which allow you to run batches of parallel code to drastically improve simulation time. You will first have to get setup with a CARC account, which can be down through your PI.

1. Open terminal and SSH using your credentials. This allows you to remotely work on the USC computers 

`ssh <username>@discovery.usc.edu`

Enter password

2. Set up a conda environment by running the following commands in terminal

`module purge`
`module load conda`

`mamba init bash`
`source ~/.bashrc`

`mamba create --prefix /project/mel_333/<carc-env-name>`

3. Activate conda environment

`mamba activate /project/mel_333/<carc-env-name>`

4. Install packages

`mamba install neuron, matplotlib, numpy, etc.`

You may also be able to use pip here. You'll want to download all the same packages (including neuron and jupyter). You shouldn't have to download XQuartz or command-line tools, though.

5. Install Jupyter kernel

`mamba install jupyter`
`mamber install ipykernel`

6. Open Jupyter Lab through CARC On Demand

Visit this [link](https://hpcaccount.usc.edu)

## Submitting Batch Jobs

1. Create a .sh file.

For batch jobs, you will want to create a .sh file which specifies how the system will run your Python scripts. To do so, you will want to create a new file using **pico**, which is a simple command-line text editor

In terminal, run

`pico <batch-file-name>.sh`

This will open pico where you will be able to edit

2. Use the following format

```
#!/bin/bash

#SBATCH --account=mel_333
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=0:30:00            # Specifies the max amoutn of time the job will running before terminating

module purge
eval "$(conda shell.bash hook)"
conda activate /project/mel_333/<env-name>   # Specifies environment path
module load python/3.9.12
cd /project/mel_333/<carc-directory>         # Specifies path for where your python file is stored      

python <python-file>.py
```

Note: You can use this exact format (i.e you shouldn't have to change much above). You just have to change the <env-name> to be the same path and name as your environment that you setup in step one. You will also have to change the cd command in to be the right file path to your python file, which should be called in the last line

3. Use ^+ O to save the .sh file
4. This will save to the Home Directory
5. Use the following to execute batch

`sbatch batch-name.sh`

You can track the progress using the CARC OnDemand link above.
