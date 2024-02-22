
local_curves = {

    'Iberia': [40.4, -3.7, 'psvc_iberia_Molina-Cardin et al_2018.dat'],

    'Iberia (Iron Age)': [40.4, -3.7, 'iberia_osete_2020.txt'],

    'France (directional)': [48.9, 2.3, 'France_LeGoff2020.dat'],

    'France (intensity)': [47.03, 4.83, 'France_LeGoff2020.dat'],

    'Western Europe (directional 1st millennium BC, Hervé et al. 2013)': [48.9, 2.3, 'we_dir1kbc.txt'],

    'Western Europe (intensity 1st millennium BC, Hervé et al. 2013)': [48.9, 2.3, 'we_int1kbc.txt'],

    'Italy': [42.45, 12.03, 'Italy_TemaLanos2020.dat'],

    'Great Britain (directional)': [52.44, -1.65, 'UK_Batt2017.dat'],

    'Bulgaria': [42.7, 23.3, 'bulgaria2014.dat'],

    'Cyprus': [35.17, 33.36, 'chipre.dat'],

    'Azores': [38.5, -28, 'azores.dat'],

    'Europe Neolithic (directional)': [43.0, 11.0, 'Directional Neolithic.dat'],

    'Hawaii': [20.24, -155.87, 'Hawaii_Tema2017.dat'], #OJO, no pone las coordenadas en la ref.

    'New Zealand': [-40, 175, 'NZ11k.dat'],

    'New Zealand (high resolution last 1 kyr)': [-40, 175, 'NZ1k.dat'],

    'New One':[],

}


global_models = {

    'SHA.DIF.14k': ['shadif14k_c.dat', 'shadif14k_ec.dat', 'shadif14k_t.dat'],
    'SHAWQ2k': ['shawq2k_c.dat', 'shawq2k_ec.dat', 'shawq2k_t.dat'],
    'SHAWQ-IronAge': ['shawqIA_c.dat', 'shawqIA_ec.dat', 'shawqIA_t.dat'],
    'ArchKalmag14k': ['archkalmag14k_c.dat', 'archkalmag14k_ec.dat', 'archkalmag14k_t.dat'],
    'BIGMUDI4k.1': ['bigmudi4k1_c.dat', 'bigmudi4k1_ec.dat', 'bigmudi4k1_t.dat'],
    'arch3k1': ['arch3k1_c.dat', 'arch3k1_ec.dat', 'arch3k1_t.dat'],
    'cals3k4': ['cals3k4_c.dat', 'cals3k4_ec.dat', 'cals3k4_t.dat'],
    'cals10k1b': ['cals10k1b_c.dat', 'cals10k1b_ec.dat', 'cals10k1b_t.dat'],

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
