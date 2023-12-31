# ArchaeoPyDating (under construction)
ArchaeoPyDating is a Python tool designed for archaeomagnetic dating. This tool is an updated release of the Matlab tool archaeo_dating (Pavón-Carrasco et al., 2011). ArchaeoPyDating is available in two formats:
- **Online Form:** For a quick and user-friendly experience, access our tool at http://pc213fis.fis.ucm.es:8080. Simply input your archaeomagnetic data, select a reference curve and a confidence limit, and click submit. The tool will process your data online and present your dating results within seconds. This format requires no coding knowledge or installations.
- **Python Module:** For users with programming skills or those seeking a more in-depth and flexible analysis. You can import this module into your Python scripts for more complex studies, automation, and integration with other libraries. It’s also compatible with notebooks for comprehensive coding studies.


##  Citation ##

## Installation ##

You can install the ArchaeoPyDating module by directly running the following command in your terminal (please ensure you have Python, pip, and git installed on your system):
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

- [Original paper](http:doi) (when available)

- [ejemplo.pdf](ejemplo.pdf): User's manual (when available)

- [APD_example.ipynb](APD_example.ipynb): A Jupyter notebook showing examples of the module's use

## Contact ##

Mario Serrano Sánchez-Bravo  
Dpto. Física de la Tierra y Astrofísica. Facultad de Ciencias Físicas.  
Universidad Complutense de Madrid (UCM)  
marioser@ucm.es  

