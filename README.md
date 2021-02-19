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

### macOS

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
