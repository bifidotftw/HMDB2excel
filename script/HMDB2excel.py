import pandas as pd
import xml.etree.ElementTree as ET



path = '../data/'


def iter_metabolites(filename):
    ctx = ET.iterparse(filename)
    for _, metabolite in ctx:
        if metabolite.tag != '{http://www.hmdb.ca}metabolite':
            continue
        yield metabolite
        metabolite.clear()


# Initialize dataframes for later appending
df_HMDB_name = pd.DataFrame()
df_HMDB = pd.DataFrame()
df_HMDB_mean = pd.DataFrame()


# For every metabolite in XML
for metabolite in iter_metabolites(path + 'HMDB_metabolites_v4.0.xml'):
    # Reset df for new iteration
    df_metabolite = pd.DataFrame()

    HMDB_id = metabolite.find('{http://www.hmdb.ca}accession').text
    name = metabolite.find('{http://www.hmdb.ca}name').text
    df_HMDB_name_x = pd.DataFrame({'HMDB_id': [HMDB_id], 'compound': [name]})
    df_HMDB_name = df_HMDB_name.append(df_HMDB_name_x)
    print(HMDB_id + ': ' + name)

    # For every entry in 'normal concentrations'
    for concentration in metabolite.find('{http://www.hmdb.ca}normal_concentrations'):
        conc = concentration.find('{http://www.hmdb.ca}concentration_value').text
        # For some reason some entries have no concentration, skip those
        if conc != None:
            '''
            HMDB entries are a mess
            If there is no space it might just be a single value
            If anything other than numbers is present it will be set None in the next step
            '''
            try:
                avg_value = conc.split(' ')[0]
            except:
                avg_value = conc
            '''
            Some concentrations are in the form of min-max
            These are skipped because we don't know where the mean is
            '''
            try:
                avg_value = float(avg_value)
            except:
                # Before adding to df this variable is checked and values are not added if avg_value == None
                avg_value = None
            '''
            Deviations are not in consistend format
            Values are sorted in value (+- x) and range (x-y)
            Rest is discarded

            Try because there might be no deviation specified
            '''
            try:
                if '+' in conc:
                    dev_value = conc.split(' ')[2]
                    dev_value = float(dev_value)
                    deviation = 'value'
                elif '(' in conc:
                    dev_value = conc.split(' ')[1]
                    dev_value = dev_value[1:-1]
                    deviation = 'range'
                else:
                    dev_value = None
                    deviation = None
            except:
                dev_value = None
                deviation = None
            '''
            Try and except in case information is missing
            '''
            try:
                unit = concentration.find('{http://www.hmdb.ca}concentration_units').text
            except:
                unit = None
                print(HMDB_id + ': no unit. Skipping!')
            try:
                sample = concentration.find('{http://www.hmdb.ca}biospecimen').text
            except:
                sample = None
            try:
                age = concentration.find('{http://www.hmdb.ca}subject_age').text
            except:
                age = None
            try: 
                sex = concentration.find('{http://www.hmdb.ca}subject_sex').text
            except:
                sex = None
            try:
                comment = concentration.find('{http://www.hmdb.ca}comment').text
            except:
                comment = None
            '''
            If there is a value, a unit and it was measured in blood convert it to df
            Else do nothing
            '''
            if  avg_value != None and unit != None and sample == 'Blood':
                df_measurement = pd.DataFrame({'HMDB_id': [HMDB_id], 'metabolite': [name], 'avg_value': [avg_value], 'dev_value': [dev_value],  'deviation': deviation, 'unit': [unit], 'sample': [sample], 'age': [age], 'sex': [sex], 'comment': [comment]})
                '''
                Collect all entries of a metabolite in single df
                '''
                df_metabolite = df_metabolite.append(df_measurement, ignore_index=True)
            else:
                pass
    # Collect all entries of every metabolite in single df
    df_HMDB = df_HMDB.append(df_metabolite)

    # Do statistics for every metabolite and collect in single dataframe
    if df_metabolite.empty == False:
        mean = df_metabolite['avg_value'].mean()
        SEM = df_metabolite['avg_value'].sem()
        n = df_metabolite.shape[0]
        metabolite_mean = pd.DataFrame({'HMDB_id': [HMDB_id], 'metabolite': [name], 'mean': [mean], 'SEM': [SEM], 'n': n})
        df_HMDB_mean = df_HMDB_mean.append(metabolite_mean, ignore_index=True)

df_HMDB_name.to_excel(path + 'HMDB_id2name.xlsx', index=False)
df_HMDB.to_excel(path + 'HMDB_blood.xlsx', index=False)
df_HMDB_mean.to_excel(path + 'HMDB_blood_avg.xlsx', index=False)
