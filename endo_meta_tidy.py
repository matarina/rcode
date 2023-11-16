
#tidy  ucec and tissue micro array meta data for clam training
import pandas as pd
import os

# |%%--%%| <vcASrZrGiG|HMf64IVk6A>

df = pd.read_csv('/data/dk/ucec/CLAM/dataset_csv/clinical.tsv', sep='\t')
cli = pd.read_csv('/data/dk/ucec/CLAM/dataset_csv/gdc_manifest.2023-09-01.txt',sep='\t')

# |%%--%%| <HMf64IVk6A|iF7yAsO80Y>

cli['case_submitter_id'] = cli['filename'].str.extract(r'(^[A-Z0-9]{4}-[A-Z0-9]{2}-[A-Z0-9]{4})')
df = df.drop_duplicates(subset='case_submitter_id',keep='first')
df = df.merge(cli, on='case_submitter_id', how='left')

# |%%--%%| <iF7yAsO80Y|ugqgmQroHb>

replacements = {
    r'Stage IV.*': "Stage_IV",
    r'Stage III.*': "Stage_III",
    r'Stage II.*': "Stage_II",
    r'Stage I.*': "Stage_I"
}
for pattern ,replace in replacements.items():
    df['figo_stage'] = df['figo_stage'].str.replace(pattern,replace)
df = df.drop('case_id',axis=1)

# |%%--%%| <ugqgmQroHb|XSDTUD7WbZ>

df['slide_id'] = df['filename'].str.replace('.svs','')
df.rename(columns={'case_submitter_id':'case_id','figo_stage':'label'},inplace=True)

# |%%--%%| <XSDTUD7WbZ|EubfO3iZ9I>

df[['label','case_id','slide_id']].to_csv('/data/dk/ucec/CLAM/dataset_csv/ucec_meta.csv', index=False)

# |%%--%%| <EubfO3iZ9I|Qf3OG0vaSB>

### tma meta rehsape below
endo = pd.read_csv('/data/dk/endometrium/endo_tma_meta.csv',sep = ',')

# |%%--%%| <Qf3OG0vaSB|hIT68Rtfle>

endo['slide_id'] = endo['序号'].apply(lambda x: x[0] + x[1:].lstrip('0'))
replacements = {
    r'T1.*': "Stage_I",
    r'T4.*': "Stage_IV",
    r'T3.*': "Stage_III"
}
endo['label'] = endo['TNM']
for pattern, replace in replacements.items():
    print(pattern)
    endo['label'] = endo['label'].str.replace(pattern,replace)
endo = endo.rename(columns={'病理号':'case_id'})
endo = endo[['case_id','slide_id','label']]


# |%%--%%| <hIT68Rtfle|51k5K7Rwil>

endo.to_csv('/data/dk/ucec/CLAM/dataset_csv/endo_meta.csv',index=False)

# |%%--%%| <51k5K7Rwil|RxnnXrA5o1>


