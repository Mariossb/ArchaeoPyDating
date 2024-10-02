
local_curves = {

    'Iberia': [40.4, -3.7, 'psvc_iberia_Molina-Cardin et al_2018.dat'],

    'Iberia (Iron Age)': [40.4, -3.7, 'iberia_osete_2020.txt'],

    'France (directional)': [48.9, 2.3, 'France_LeGoff2020.dat'],

    'France (intensity)': [47.03, 4.83, 'FranceIntensityGenevey2021_F6b.txt'],

    'Western Europe (directional 1st millennium BC, Hervé et al. 2013)': [48.9, 2.3, 'we_dir1kbc.txt'],

    'Western Europe (intensity 1st millennium BC, Hervé et al. 2013)': [48.9, 2.3, 'we_int1kbc.txt'],

    'Italy': [42.45, 12.03, 'Italy_TemaLanos2020.dat'],

    'Great Britain (directional)': [52.44, -1.65, 'UK_Batt2017.dat'],

    'Bulgaria': [42.7, 23.3, 'bulgaria2014.dat'],

    'Cyprus': [35.17, 33.36, 'chipre.dat'],

    'Azores': [38.5, -28, 'azores.dat'],

    'Europe Neolithic (directional)': [43.0, 11.0, 'Directional Neolithic.dat'],

    'Hawaii': [20.24, -155.87, 'Hawaii_Tema2017.dat'],

    'New Zealand': [-40, 175, 'NZ11k.dat'],

    'New Zealand (high resolution last 1 kyr)': [-40, 175, 'NZ1k.dat'],

    'Mexico (directional)': [19.4, -99.1, 'Mexico_dir.dat'],

    'Mexico (intensity)': [19.4, -99.1, 'Mexico_int.txt'],

    'LAC.v.2.0 (Levant intensity)': [31.78, 35.21, 'LAC_V_2_0.txt'],

    'New One':[],

}


global_models = {

    'SHA.DIF.14k': ['shadif14k_c.dat', 'shadif14k_ec.dat', 'shadif14k_t.dat'],
    'SHAWQ2k': ['shawq2k_c.dat', 'shawq2k_ec.dat', 'shawq2k_t.dat'],
    'SHAWQ-IronAge': ['shawqIA_c.dat', 'shawqIA_ec.dat', 'shawqIA_t.dat'],
    'ArchKalmag14k.r': ['archkalmag14k_c.dat', 'archkalmag14k_ec.dat', 'archkalmag14k_t.dat'],
    'BIGMUDI4k.1': ['bigmudi4k1_c.dat', 'bigmudi4k1_ec.dat', 'bigmudi4k1_t.dat'],

}


regional_models = {

    'SCHA.DIF.4k': 'scha_dif_4k_archaeo_dating.dat',
    'SCHAFRICA.DIF.4k': 'schafrica_dif_4k_archaeo_dating.dat',

}


curves = {

'Select one': [''],

'Local Curve': list(local.keys()),

'Regional Model': list(regional_models.keys()),

'Global Model': list(global_models.keys()),


}


curves_references = {
    "Local PSVCs": {
        'Iberia': 'Molina-Cardín et al. (2018) [D,I,F] [1000 BC - 1900 AD]',
        'Iberia (Iron Age)': 'Osete et al. (2020) [D,I,F] [1100 BC - 100 BC]',
        'France (directional)': 'Le Goff et al. (2020) [D,I] [70 BC - 1700 AD]',
        'France (intensity)': 'Genevey et al. (2021) [F] [300 AD - 1850 AD]',
        'Western Europe (directional 1st millennium BC, Hervé et al. 2013)': 'Hervé et al. (2013a) [D,I] [1500 BC - 200 AD]',
        'Western Europe (intensity 1st millennium BC, Hervé et al. 2013)': 'Hervé et al. (2013b) [F] [1500 BC - 200 AD]',
        'Italy': 'Tema and Lanos (2021) [D,I,F] [900 BC - 1980 AD] (F only for [700 BC - 1800 AD])',
        'Great Britain (directional)': 'Batt et al. (2017) [D,I] [6000 BC - 1980 AD]',
        'Bulgaria': 'Kovacheva et al. (2014) [D,I,F] [6000 BC - 1800 AC]',
        'Cyprus': 'Tema et al. (2021) [D,I] [2000 BC - 1900 AC]',
        'Azores': 'Béguin et al. (2020) [D,I,F] [500 BC - 2000 AD]',
        'Europe Neolithic (directional)': 'Carrancho et al. (2013) [D,I] [6000 BC - 1000 AD]',
        'Hawaii': 'Tema et al. (2017) [D,I,F] [8000 BC - 1950 AD]',
        'New Zealand': 'Turner and Corkill (2023) [D,I,F] [9325 BC - 1875 AD]',
        'New Zealand (high resolution last 1 kyr)': 'Turner and Corkill (2023) [D,I,F] [1025 AD - 2000 AD]',
        'Mexico (directional)': 'García-Ruiz et al. (2022)[D,I] [1472 BC - 1947 AD]',
        'Mexico (intensity)': 'García-Ruiz et al. (2021)[F] [1600 BC - 1946 AD]',
        'LAC.v.2.0 (Levant intensity)': 'Hassul et al. (2024) [F] [3200 BC - 550 AD]',
    },
    "Global Models": {
        'SHA.DIF.14k': 'SHA.DIF.14k (Pavón-Carrasco et al., 2014) [12000 BC - 1900 AD]',
        'SHAWQ2k': 'SHAWQ2k (Campuzano et al., 2019) [100 BC - 1900 AD]',
        'SHAWQ-IronAge': 'SHAWQ-IronAge (Osete et al., 2020) [1300 BC - 0 AD]',
        'ArchKalmag14k.r': 'ArchKalmag14k.r (Schanner et al., 2022) [12000 BC - 1950 AD]',
        'BIGMUDI4k.1': 'BIGMUDI4k.1 (Arneitz et al., 2019) [2000 BC - 2000 AD]',
    },
    "Regional Models": {
        'SCHA.DIF.4k': 'SCHA.DIF.4k (Pavón-Carrasco et al., 2021) [2000 BC - 1900 AD]',
        'SCHAFRICA.DIF.4k': 'SCHAFRICA.DIF.4k (Di Chiara and Pavón-Carrasco, 2022) [2000 BC - 1900 AD]',
    }
}

# Función para imprimir todas las curvas unificadas
def available_psvc():
    for category, curves in curves_references.items():
        print(f"\n{category}:")
        for curve, reference in curves.items():
            print(f"- {curve}: {reference}")


