# Reg-Seq Returns - The Whole Genome

## Overview
This repository is the beginning the efforts to improve the original [Reg-Seq](https://github.com/RPGroup-PBoC/RegSeq) method to usable on a genome wide scale. This repository contains a software module, which can be found in the `software_module` folder. To install it, navigate into the folder, and install required packages by running `pip install -r requirements.txt`. Then, simply install the module by running `pip install -e .`, which installs the package in the editable mode, so modifications can be made without reinstalling it.


## Layout 

The repository is split into seven main directories, many of which have subdirectories. This structure has been designed to be easily navigable by humans and computers alike, allowing for rapid location of specific files and instructions. Within each directory is a `README.md` file which summarizes the purpose of that directory as well as some examples where necessary. This structure may not be perfect for your intended us and may need to be modified. Each section is briefly described below.

### **`code`**
Where all of the *executed* code lives. This includes pipelines, scripts, and figure files.
 * **`processing`**: Any code used to *transform* the data into another type should live here. This can include everything from parsing of text data, image segmentation/filtering, or simulations.
 * **`analysis`**: Any code to *draw conclusions* from an experiment or data set. This may include regression, dimensionality reduction, or calculation of various quantities.
 * **`exploratory`**: A sandbox where you keep a record of your different approaches to transformation, interpretation, cleaning, or generation of data.
 * **`figures`**: Any code used to generate figures for your finished work, presentations, or for any other use.

### **`data`**
All raw data collected from your experiments as well as copies of the transformed data from your processing code.

### **`miscellaneous`**
Files that may not be code, but are important for reproducibility of your findings.
* **`protocols`**: A well annotated and general description of your experiments. These protocols should be descriptive enough for someone to follow your experiments independently.
* **`materials`**: Information regarding the materials used in your experiments or data generation. This could include manufacturer information, records of purity, and/or lot and catalog numbers.
* **`software details`**: Information about your computational environment that are necessary for others to execute your code. This includes details about your operating system, software version and required packages.

### **`software_module`**
Contains a python module. Installed by navigating into this folder in the command line and using `pip install -e .`. Necessary to run most of the code used in the `code` folder. 

### **`templates`**
Files that serve as blank templates that document the procedures taken for each experiment, simulation, or analysis routine.



# License Information


# License Information
<img src="https://licensebuttons.net/l/by/3.0/88x31.png"> This work is
licensed under a [Creative Commons CC-BY 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/). All
software is issued under the standard MIT license which is as follows:

```
Copyright 2020, The authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
