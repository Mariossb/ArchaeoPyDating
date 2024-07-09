# ArchaeoPyDating
ArchaeoPyDating (Serrano et al., 2024) is a Python tool designed for archaeomagnetic dating. This tool is an updated release of the Matlab tool archaeo_dating (Pavón-Carrasco et al., 2011). ArchaeoPyDating is available in two formats:
- **Online Form:** For a quick and user-friendly experience, access our tool at http://pc213fis.fis.ucm.es:8080. Simply input your archaeomagnetic data, select a reference curve and a confidence limit, and click submit. The tool will process your data online and present your dating results within seconds. This format requires no coding knowledge or installations.
- **Python Module:** For users with programming skills or those seeking a more in-depth and flexible analysis. You can import this module into your Python scripts for more complex studies, automation, and integration with other libraries. It’s also compatible with notebooks for comprehensive coding studies.


##  How to cite ##

Serrano, M., Pavón-Carrasco, F.J., Campuzano, S.A., & Osete, M.L. (2024). ArchaeoPyDating: A new user-friendly release for archaeomagnetic dating. Archaeometry, 1–14. https://doi.org/10.1111/arcm.13009

## Using the tool in Google Colab ##

You can use ArchaeoPyDating directly in a Google Colab notebook. Here's how:

1. **Open a new book**

2. **Install the module in the current notebook:**
   ```
   ! pip install -q git+https://github.com/Mariossb/ArchaeoPyDating.git
   ```

4. **Import the module:**
      ```python
   import ArchaeoPyDating.ArchaeoPyDating as apd
   ```
5. **Usage**:  Now, you can utilize the objects and functions provided by the ArchaeoPyDating module within this Google Colab notebook.

## Installation ##

You can install the ArchaeoPyDating module in your own system by directly running the following command in your terminal (please ensure you have Python, pip, and git installed on your system):
```
pip install git+https://github.com/Mariossb/ArchaeoPyDating.git
```
If you have a virtual environment activated, the module will be installed in that virtual environment. If you do not have a virtual environment activated, it will be installed in the global Python environment on your system.

If you want to create a directory within which you can create a virtual Python environment for your archaeomagnetic dating projects, follow these steps:

1. **Create a new directory**: Open your terminal, create a new directory, and/or navigate to it:
   ```
   mkdir mydirectoryname
   cd mydirectoryname
   ```

2. **Create a virtual Python environment**: 
    ```
    python -m venv myenvname
    ```

3. **Activate the virtual environment**: 
   ```
   myenvname\Scripts\activate  # On Windows systems
   ```   
   ```
   source myenvname/bin/activate  # On Unix/Linux-based systems
   ```

4. **Intall the module along with all dependencies**: 
   ```
   pip install git+https://github.com/Mariossb/ArchaeoPyDating.git
   ```

5. **Use it**:  Now, you can use the objects and functions of the ArchaeoPyDating module within this virtual environment. To import them into your .py files, include the following line at the beginning (see this [notebook](APD_example.ipynb) for examples):
   ```python
   import ArchaeoPyDating as apd
   ```



##  Documentation ##

- [Original paper](https://doi.org/10.1111/arcm.13009)

- [Notebook showing examples of the Python module's use](APD_example.ipynb)

## Contact ##

Mario Serrano Sánchez-Bravo  
Dpto. Física de la Tierra y Astrofísica. Facultad de Ciencias Físicas.  
Universidad Complutense de Madrid (UCM)  
marioser@ucm.es  || paleomag@ucm.es

