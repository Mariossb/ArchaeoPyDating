# ArchaeoPyDating
ArchaeoPyDating is a Python tool designed for archaeomagnetic dating. This tool is an updated release of the Matlab tool archaeo_dating (Pavón-Carrasco et al., 2011). ArchaeoPyDating is available in two formats:
- **Online Form:** For a quick and user-friendly experience, access our tool at http://pc213fis.fis.ucm.es:8080. Simply input your archaeomagnetic data, select a reference curve and a confidence limit, and click submit. The tool will process your data online and present your dating results within seconds. This format requires no coding knowledge or installations.
- **Python Module:** For users with programming skills or those seeking a more in-depth and flexible analysis. You can import this module into your Python scripts for more complex studies, automation, and integration with other libraries. It’s also compatible with notebooks for comprehensive coding studies.


##  Instructions ##

- Data: This class represents the undated archaeological artefact or lava flow. It
accepts the following key arguments: declination (D=), inclination (I=), α95 angle
(a95=), intensity (F=), intensity standard deviation (eF=), site latitude (lat=),
longitude (lon=) and sitename (sitename=). Angles are assumed to be in decimal
degrees and intensities in micro-Teslas. The class includes three methods:

   - Data.cvp(nlat, nlon): Translates the directional data to a new lati-
   tude/longitude location (nlat/nlon) using the conversion via pole (CVP)
   method.
  
   - Data.vadm(nlat): Translates the intensity data to a new latitude (nlat)
   using the virtual axial dipole moment (VADM) method.
  
   - Data.vdm: Translates the intensity data using the virtual dipole moment
   (VDM) method. This method should be used after applying the CVP method.

- Curve: This class represents the paleosecular variation curve (PSVC) used in the
dating process. Its construction dependas on the desired type of curve. There are
three options:

   - Curve(local = input1): Constructs a classical PSVC. Available options for
   input1 are: ’Iberia’, ’France’, ...

   - Curve(gmodel = input1, lat = input2, lon = input3): Constructs
   a synthetic PSVC derived from a global paleomagnetic model. Available
   options for input1 are: ’SHA.DIF.14k’, ’SHAWQ2k’, ’SHAWQ-IronAge’, ’ArchKalmag14k’, ’BIGMUDI4k.1’, ’arch3k1’, ’cals3k4’, ’cals10k1b’. Latitude
   (input2 ) and longitude (input3 ) where the PSVC needs to be synthesized are
   required.

   - Curve(rmodel = input1, lat = input2, lon = input3): Constructs a
   PSVC derived from a regional paleomagnetic model. Available options for
   input1 are: ’SCHA.DIF.4k’, ’SCHAFRICA.DIF.4k’. Latitude (input2 ) and
   longitude (input3 ) where the PSVC needs to be synthesized are also required.

   The Curve class has one method:

   - Curve.int temp(tmin, tmax): Temporally delimits the curve by specifying
   a minimum and maximum time (tmin/tmax).

- Dating: This class combines the two previous classes to perform archaeomagnetic
dating. It is constructed as Dating(Data, Curve) and includes five methods:

   - Dating.datingX(element). Returns a numpy.array with the probability
   density function (PDF) of the corresponding geomagnetic element (elements
   can be ’D’, ’I’ or ’F’ ). The time-axis of the PDF ranges from the minimum
   time of the Curve to its maximum in steps of 1 year.

   - Dating.zcomb(*arg): Returns a combined PDF by taking multiple PDFs
   (obtained from the above method) as arguments.

   - Dating.pb(PDF, p): Returns the threshold in the PDF for a given probability (p, e.g., 95%).

   - Dating.pb h(PDF, h): Returns the final date given the PDF and the
   threshold, i.e. with what is obtained in the two previous methods.

   - Dating.plot(zD= , zI= , zF= , z= , hD= , hI= , hF= , h= ). Returns
   the final figure representing the dating. The keys arguments refer to the PDFs
   in the case of z (e.g., zD is the PDF corresponding to declination, and z is the
   PDF corresponding to the total PDF), and to the thresholds in the case of h


##  Citation ##

## Installation ##

Follow these steps to install and use the ArchaeoPyDating module (please, ensure you have Python and pip installed on your system):

1. **Clone the repository**: Open your terminal, create a new directory or navigate to a desired one, and run the following command:
   ```
    git clone https://github.com/Mariossb/ArchaeoPyDating.git
   ```

2. **Navigate to the project directory**: Change your current directory to the cloned repository's directory:
    ```
    cd ArchaeoPyDating
    ```

3. **Install the required libraries**: Use pip in your terminal:
   ```
   pip install -r requirements.txt
   ```

4. **Import the module**: Create your own .py file and import the module by writting at the beggining:
   ```python
   import ArchaeoPyDating as apd
   ```

5. **Use it**: Now, you can use the objects and functions of the ArchaeoPyDating module. See this [Jupyter Notebook](APD_example.ipynb) for examples.




##  Documentation ##

- [Original paper](http:doi)

- [ejemplo.pdf](ejemplo.pdf): User's manual

- [APD_example.ipynb](APD_example.ipynb): A Jupyter notebook showing examples of the module's use

## Contact ##

Mario Serrano Sánchez-Bravo  
Dpto. Física de la Tierra y Astrofísica. Facultad de Ciencias Físicas.  
Universidad Complutense de Madrid (UCM)  
marioser@ucm.es  

