# Lecture Supplementary Material 2021

This is a collection of codes and Jupyter notebooks to supplement the 2021 ISR summer school lectures.

## Using in the Cloud
You can interact with the Jupyter notebooks in this repository using the free server, Binder. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/isrsummerschool/lecture_supplements_2021/main)


## Using Locally

If you want to use this locally, here are the basics steps:
1) clone this repository
2) set up a python environment
3) run a Jupyter Lab server

Below, we have some basic instructions for each operating system to guide you:

### Windows
#### Recommended Python and Git Installation
We recommend using a package manager to obtain several common developer tools, including git and miniconda. We'll use git to clone git repositories and miniconda to provide a Python installation. Miniconda is a free minimal installation of conda. Many Python packages depend on other system libraries and conda provides both the Python packages and the system libraries in a platform independent way. Basically this helps make using Python on Windows and macOS much easier.
1.	Download and install chocolatey. It’s a package manager for helping you install and manage other things: https://chocolatey.org/
2.	In a powershell terminal with administrative privledges, install git, bash, miniconda, and the mingw32 compiler:
    1. choco install git
    2. choco install miniconda3 --params="'/AddToPath:1'" --params="'/RegisterPython:1'"

For more information, see the documentation on chocolatey for [git](https://community.chocolatey.org/packages/git), [miniconda3](https://community.chocolatey.org/packages/miniconda3), and [anaconda3](https://community.chocolatey.org/packages/anaconda3)

##### Using Anaconda instead of miniconda
If you would prefer to install the full Anaconda instead of miniconda, replace step b. with:
b. choco install anaconda3 --params="'/AddToPath:1'"

#### Installing Python and Git Without Chocolatey
If you would prefer to install git and miniconda without using chocolatey, you can do so by following these instructions:
1. Download and install git-bash from [here](https://gitforwindows.org/)
2. Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/products/individual-b)
### macOS

**1)** Install Anaconda for mac from https://www.anaconda.com/products/individual

**2)** Set up a python environment

    Open a terminal
    > conda create -n isrschool anaconda  (or use a different name than isrschool)
    > conda activate isrschool
    
**3)** Install git

    > conda install git 
    
**4)** Clone this repository:

    $ git clone https://github.com/isrsummerschool/lecture_supplements_2021.git
    
**5)** Install the python packages that the lecture supplement notebooks need:

    $ cd lecture_supplements_2021
    $ pip install -r requirements.txt


**6)** Now you can start up a Jupyter Lab server and work with the notebooks:

    $ jupyter lab


### Linux

We only provide example commands here for Ubuntu, but the process is similar for other Linux distributions.

**1)** Install ``git`` with your package manager and clone the repository:

Ubuntu:

    $ sudo apt-get install git-all

**2)** Clone this repository:

    $ git clone https://github.com/isrsummerschool/lecture_supplements_2021.git
    
**3)** Set up a python environment (for more details read this: https://virtualenv.pypa.io/en/latest/user_guide.html)

Install the Python ``pip`` and ``virtualenv`` packages:

    $ sudo apt install python3-pip
    $ pip install virtualenv
    
and then create a virtual environment:

    $ virtualenv isrschool
    $ source isrschool/bin/activate
    
and finally, install the python packages that the lecture supplement notebooks need:

    $ cd lecture_supplements_2021
    $ pip install -r requirements.txt


**4)** Now you can start up a Jupyter Lab server and work with the notebooks:

    $ jupyter lab
