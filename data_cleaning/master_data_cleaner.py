#! /usr/bin/env python3
import numpy as np
import pandas as pd


def stage_changer(stage_string):
    """
    Input: A string specifying the stage of the given tumor.
    Output: An integer representing the stage of the given tumor.
    """
    if (stage_string == 'stage i') or (stage_string == 'stage ia') or (stage_string == 'stage ib') or (stage_string == 'stage ic') or (stage_string == 'i') or (stage_string == 'ia') or (stage_string == 'ib') or (stage_string == 'ic'):
        return 1
    elif (stage_string == 'stage ii') or (stage_string == 'stage iia') or (stage_string == 'stage iib') or (stage_string == 'stage iic') or (stage_string == 'ii') or (stage_string == 'iia') or (stage_string == 'iib') or (stage_string == 'iic'):
        return 2
    elif (stage_string == 'stage iii') or (stage_string == 'stage iiia') or (stage_string == 'stage iiib') or (stage_string == 'stage iiic') or (stage_string == 'iii') or (stage_string == 'iiia') or (stage_string == 'iiib') or (stage_string == 'iiic'):
        return 3
    elif (stage_string == 'stage iv') or (stage_string == 'stage iva') or (stage_string == 'stage ivb') or (stage_string == 'stage ivc') or (stage_string == 'iv') or (stage_string == 'iva') or (stage_string == 'ivb') or (stage_string == 'ivc'):
        return 4
    elif (stage_string == 'stage x') or (stage_string == 'x'):
        return 10


def death_finder(row, death_columns, follow_up_columns):
    """
    Input: A row of data representing a patient record, a list of columns to
    check for days to death, and a list of columns to check for days to last
    followup.
    Output: A tuple of integers. The first integer being a 0 if the patient is
    recorded as dead at any followup and a one otherwise. The second integer
    being the day of the last followup (whether the patient was alive or dead).
    The third integer being the age of the patient at the last followup
    (whether the patient was alive or dead).
    """
    # Default value patient is alive. Search the entries to see if the patient
    # is recorded as dead at any point.
    aod = 1
    for column in death_columns:
        if float(row[column]) > 0:
            aod = 0
    # If the patient is recorded as dead, use last days to death entry as being
    # the last followup. If the patient is alive, use the last days to followup
    # entry as the day of the last followup.
    if aod == 0:
        death_follow_ups = [float(row[column]) for column in death_columns]
        followup = max(death_follow_ups)
    elif aod == 1:
        alive_follow_ups = [float(row[column]) for column in follow_up_columns]
        followup = max(alive_follow_ups)
    age = followup - float(row['patient.days_to_birth'])
    return aod, followup, age


def cleaner(clin_data, gene_data, stage_column, death_columns, follow_up_columns,
            clean_clin_file, clean_gene_file):
    """
    Input: String specifying the text file of the clinical data, string
    specifying the text file of the normalized gene sequencing data, heading to
    check for the stage of the tumor, list of headings to check for patient
    death in the clinical data, list of headings to check for last follow up
    with the patient in the clinical data, string specifying the name of the
    cleaned clinical data to be written, and string specifying the name of the
    cleaned gene expression data file to be written.
    Output: CSV files with samples from neighboring tissue removed, rows of
    patient data in alphabetical order by patient barcode, only including those
    patients for which both clinical and gene expression data are available.
    Return None.
    """
    # Read the text file.
    clin_clean = pd.read_csv(clin_data, sep='\t', header=None)
    clin_clean = clin_clean.set_index(0)
    # Transpose the dataframe so that each column represents a kind of data,
    # rather than a patient as in the orignal dataset.
    clin_clean = clin_clean.transpose()
    # Alphabetize by patient barcode.
    clin_clean = clin_clean.sort_values(by='patient.bcr_patient_barcode')
    # Read the text file.
    gene_clean = pd.read_csv(gene_data, sep='\t', header=None)
    gene_clean = gene_clean.set_index(0)
    # Transpose the dataframe so that each column represents a kind of data,
    # rather than a patient as in the orignal dataset.
    gene_clean = gene_clean.transpose()
    # Alphabetize by patient barcode.
    gene_clean = gene_clean.sort_values(by='Hybridization REF')
    # Remove samples taken from neighboring tissue.
    for index, data in gene_clean.iterrows():
        if data['Hybridization REF'][13:16] == '11A':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        elif data['Hybridization REF'][13:16] == '11B':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        elif data['Hybridization REF'][13:16] == '02A':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        elif data['Hybridization REF'][13:16] == '02B':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        elif data['Hybridization REF'][13:16] == '05A':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        elif data['Hybridization REF'][13:16] == '06A':
            gene_clean = gene_clean.drop(labels=index, axis=0)
        else:
            pass
    # Remove from the gene expression data any samples for which there is no
    # matching record in the clinical data.
    for index, data in gene_clean.iterrows():
        if data['Hybridization REF'][:12].lower() in clin_clean ['patient.bcr_patient_barcode'].to_numpy():
            pass
        else:
            gene_clean = gene_clean.drop(labels=index, axis=0)
    # Create an array consisting of lower case values of the patient barcodes
    # in the gene expression data.
    hybrid_ref_pre = []
    for code in gene_clean['Hybridization REF'].to_numpy():
        hybrid_ref_pre += [code[:12]]
    # Remove from the clinical data any patients for which there is no matching
    # record in the gene expression data.
    for index, data in clin_clean.iterrows():
        if data['patient.bcr_patient_barcode'].upper() in hybrid_ref_pre:
            pass
        else:
            clin_clean = clin_clean.drop(labels=index, axis=0)
    # Compare the patient barcodes present in the gene expression data and the
    # clinical data to ensure every sample matches
    # hybrid_ref_pre = np.array(list(map(lambda x: x.lower(), hybrid_ref_pre)))
    # print(hybrid_ref_pre == clin_clean['patient.bcr_patient_barcode'].to_numpy())
    # Create new columns in the dataframe containing the stage of the tumor
    # as an integer whether the patient was alive or dead at
    # last contact, the days to last contact, and the patient's age at last
    # contact.
    # clin_clean['Integer Stage'] = list(map(lambda x: stage_changer(x), clin_clean[stage_column])) #Comment out if stage not present
    clin_clean['Alive or Dead'] = clin_clean.apply(lambda x: death_finder(x, death_columns, follow_up_columns)[0], axis=1)
    clin_clean['Last Contact'] = clin_clean.apply(lambda x: death_finder(x, death_columns, follow_up_columns)[1], axis=1)
    clin_clean['Age at Last Contact'] = clin_clean.apply(lambda x: death_finder(x, death_columns, follow_up_columns)[2], axis=1)
    # Transpose the dataframes so the data will be in the same format as in the
    # original text files.
    clin_clean = clin_clean.transpose()
    gene_clean = gene_clean.transpose()
    # Write the cleaned data to csv files.
    clin_clean.to_csv(clean_clin_file, header=False)
    gene_clean.to_csv(clean_gene_file, header=False)
    return

# Run the cleaner on each tumor type here for ease of data entry.
cleaner('../2022_raw_data/clinical/GBMLGG.clin.merged.txt',
       '../2022_raw_data/rna_exp/GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',
       'patient.stage_event.clinical_stage',
       ['patient.days_to_death',
        'patient.follow_ups.follow_up.days_to_death',
        'patient.follow_ups.follow_up-2.days_to_death',
        'patient.follow_ups.follow_up-3.days_to_death'],
       ['patient.days_to_last_followup',
        'patient.follow_ups.follow_up.days_to_last_followup',
        'patient.follow_ups.follow_up-2.days_to_last_followup',
        'patient.follow_ups.follow_up-3.days_to_last_followup'],
       '../2022_processed_data/GBMLGG_clin_clean.csv', '../2022_processed_data/GBMLGG_gene_clean.csv')

# Join and invert processed data on patient ID for collaborations.
import numpy as np
import pandas as pd
shared_data_genetic = pd.read_csv("2022_processed_data/GBMLGG_gene_clean.csv",header=None, index_col=0).transpose()
shared_data_clinical = pd.read_csv("2022_processed_data/GBMLGG_clin_clean.csv",header=None,index_col=0).transpose()
shared_data_complete = shared_data_genetic
shared_data_complete['Alive or Dead'] = shared_data_clinical['Alive or Dead']
shared_data_complete['Last Contact'] = shared_data_clinical['Last Contact']
shared_data_complete['Age at Last Contact'] = shared_data_clinical['Age at Last Contact']
print(shared_data_complete.head())
# complete_data = pd.concat([join_data_clinical, join_data_genetic]).transpose()
# complete_data.to_csv('2022_processed_data/GBMLGG_complete_clean.csv')