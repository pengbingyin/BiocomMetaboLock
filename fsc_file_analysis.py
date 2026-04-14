import fcsparser

file_path = r"C:\Users\pengb\PycharmProjects\FPP_FACS\data\evolution1\Strain evolution 1-1 and screening\Specimen_001_FS_004.fcs"

print("=== Reading FCS file ===")
meta, data = fcsparser.parse(file_path, reformat_meta=True)

print("\n=== FULL METADATA HEADER (all keys and values) ===")
for key in sorted(meta.keys()):
    value = meta[key]
    if isinstance(value, bytes):
        try:
            value = value.decode('utf-8', errors='ignore')
        except:
            value = str(value)
    print(f"{key:35} : {value}")

print("\n=== DATAFRAME PREVIEW ===")
print(data.head())

print("\n=== COLUMN NAMES ===")
print(list(data.columns))

print("\n=== NUMBER OF EVENTS ===")
print(len(data))