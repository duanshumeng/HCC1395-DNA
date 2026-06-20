import mysql.connector
cnn = mysql.connector.connect(user='root',passwd='Duan1996@',database='MassBank')
cursor = cnn.cursor()
query='SELECT ms_compound.compound_id,exactmass, inchi, inchikey,cas,synonym,spectrum_name,ms_level,polarity,precursor_mz_text,ionization,collision_energy_text,PK_PEAK_MZ,PK_PEAK_INTENSITY FROM ms_compound  left join synonym on ms_compound.compound_id=synonym.compound_id left join msms_spectrum on msms_spectrum.compound_id = ms_compound.compound_id left join PEAK on PEAK.RECORD = msms_spectrum.spectrum_id where ms_compound.cas = %s;'

meta_info = pd.read_excel('20221101_RefDatasets_CSA.xlsx')
for i in meta_info['CAS'].drop_duplicates().values:
    print('Processing.... ',i)
    values = [i]
    cursor.execute(query,values)
    data = cursor.fetchall()
    
    
