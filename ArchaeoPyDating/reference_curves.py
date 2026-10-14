'''Reference curves and models available in ArchaeoPyDating.

Each entry of reference_curves describes one curve or model:

- type: 'local', 'global' or 'regional'
- local curves: file (in curves/local) and lat, lon of the curve
- global models: files (coefficients, errors, epochs; in curves/gmodel)
- regional models: file (in curves/rmodel)
- reference: publication (models; local curves include it in the name)
- note: optional remark

To add a curve: put its file(s) in the right folder and add one entry.
'''

__all__ = ['reference_curves', 'local_curves', 'global_models',
           'regional_models', 'curves']


reference_curves = {

    # ---------------------------------------------------------- local curves
    'Iberia (Molina-Cardín et al. 2018)': {
        'type': 'local', 'lat': 40.4, 'lon': -3.7,
        'file': 'psvc_iberia_Molina-Cardin et al_2018.dat'},

    'Iberia (Iron Age, Osete et al. 2020)': {
        'type': 'local', 'lat': 40.4, 'lon': -3.7,
        'file': 'iberia_osete_2020.txt'},

    'France (directional, Le Goff et al. 2020)': {
        'type': 'local', 'lat': 48.9, 'lon': 2.3,
        'file': 'France_LeGoff2020.dat'},

    'France (intensity, Genevey et al. 2021)': {
        'type': 'local', 'lat': 47.03, 'lon': 4.83,
        'file': 'FranceIntensityGenevey2021_F6b.txt'},

    'Western Europe (directional 1st millennium BC, Hervé et al. 2013)': {
        'type': 'local', 'lat': 48.9, 'lon': 2.3,
        'file': 'we_dir1kbc.txt'},

    'Western Europe (intensity 1st millennium BC, Hervé et al. 2013)': {
        'type': 'local', 'lat': 48.9, 'lon': 2.3,
        'file': 'we_int1kbc.txt'},

    'Italy (Tema and Lanos 2020)': {
        'type': 'local', 'lat': 42.45, 'lon': 12.03,
        'file': 'Italy_TemaLanos2020.dat',
        'note': 'F only for 700 BC - 1800 AD'},

    'Great Britain (directional, Batt et al. 2017)': {
        'type': 'local', 'lat': 52.44, 'lon': -1.65,
        'file': 'UK_Batt2017.dat'},

    'Bulgaria (Kovacheva et al. 2014)': {
        'type': 'local', 'lat': 42.7, 'lon': 23.3,
        'file': 'bulgaria2014.dat'},

    'Cyprus (directional, Tema et al. 2021)': {
        'type': 'local', 'lat': 35.17, 'lon': 33.36,
        'file': 'chipre.dat'},

    'Azores (Béguin et al. 2020)': {
        'type': 'local', 'lat': 38.5, 'lon': -28,
        'file': 'azores.dat'},

    'Europe Neolithic (directional, Carrancho et al. 2013)': {
        'type': 'local', 'lat': 43.0, 'lon': 11.0,
        'file': 'Directional Neolithic.dat'},

    'Hawaii (Tema et al. 2017)': {
        'type': 'local', 'lat': 20.24, 'lon': -155.87,
        'file': 'Hawaii_Tema2017.dat'},

    'New Zealand (Turner and Corkill 2023)': {
        'type': 'local', 'lat': -40, 'lon': 175,
        'file': 'NZ11k.dat'},

    'New Zealand (high resolution last 1 kyr, Turner and Corkill 2023)': {
        'type': 'local', 'lat': -40, 'lon': 175,
        'file': 'NZ1k.dat'},

    'Mexico (directional, García-Ruiz et al. 2022)': {
        'type': 'local', 'lat': 19.4, 'lon': -99.1,
        'file': 'Mexico_dir.dat'},

    'Mexico (intensity, García-Ruiz et al. 2021)': {
        'type': 'local', 'lat': 19.4, 'lon': -99.1,
        'file': 'Mexico_int.txt'},

    'LAC.v.2.0 (Levant intensity, Hassul et al. 2024)': {
        'type': 'local', 'lat': 31.78, 'lon': 35.21,
        'file': 'LAC_V_2_0.txt'},

    # --------------------------------------------------------- global models
    'SHA.DIF.14k': {
        'type': 'global', 'reference': 'Pavón-Carrasco et al. 2014',
        'files': ['shadif14k_c.dat', 'shadif14k_ec.dat', 'shadif14k_t.dat']},

    'SHAWQ2k': {
        'type': 'global', 'reference': 'Campuzano et al. 2019',
        'files': ['shawq2k_c.dat', 'shawq2k_ec.dat', 'shawq2k_t.dat']},

    'SHAWQ-IronAge': {
        'type': 'global', 'reference': 'Osete et al. 2020',
        'files': ['shawqIA_c.dat', 'shawqIA_ec.dat', 'shawqIA_t.dat']},

    'ArchKalmag14k.r': {
        'type': 'global', 'reference': 'Schanner et al. 2022',
        'files': ['archkalmag14k_c.dat', 'archkalmag14k_ec.dat', 'archkalmag14k_t.dat']},

    'BIGMUDI4k.1': {
        'type': 'global', 'reference': 'Arneitz et al. 2019',
        'files': ['bigmudi4k1_c.dat', 'bigmudi4k1_ec.dat', 'bigmudi4k1_t.dat']},

    'SHAS.LIAA': {
        'type': 'global', 'reference': 'Rivera et al. 2023',
        'note': 'perturbed SHA.DIF.14k',
        'files': ['shasliaa_c.dat', 'shasliaa_ec.dat', 'shasliaa_t.dat']},

    # ------------------------------------------------------- regional models
    'SCHA.DIF.4k': {
        'type': 'regional', 'reference': 'Pavón-Carrasco et al. 2021',
        'file': 'scha_dif_4k_archaeo_dating.dat'},

    'SCHAFRICA.DIF.4k': {
        'type': 'regional', 'reference': 'Di Chiara and Pavón-Carrasco 2022',
        'file': 'schafrica_dif_4k_archaeo_dating.dat'},
}


# Views by type, in the format used by Curve, Model and RegionalModel
local_curves = {name: [c['lat'], c['lon'], c['file']]
                for name, c in reference_curves.items() if c['type'] == 'local'}
global_models = {name: c['files']
                 for name, c in reference_curves.items() if c['type'] == 'global'}
regional_models = {name: c['file']
                   for name, c in reference_curves.items() if c['type'] == 'regional'}

# Options of the online app (to be moved to the app)
local_curves['New One'] = []

curves = {
    'Select one': [''],
    'Local Curve': list(local_curves.keys()),
    'Regional Model': list(regional_models.keys()),
    'Global Model': list(global_models.keys()),
}