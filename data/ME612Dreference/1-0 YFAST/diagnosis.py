import fcsparser

file_path = r"C:\Users\pengb\PycharmProjects\FPP_FACS\data\ME612Dreference\1-0 YFAST\Specimen_001_Tube_001_001.fcs"  # ← change if needed

meta, data = fcsparser.parse(file_path, reformat_meta=True)

print("=== Type of meta ===")
print(type(meta))

print("\n=== All keys that contain 'P' and ('N' or 'V') ===")
for k in sorted(meta.keys()):
    if 'P' in str(k).upper() and ('N' in str(k).upper() or 'V' in str(k).upper()):
        print(repr(k), "→", repr(meta[k]))

print("\n=== Direct access tests ===")
print("$P4N" in meta)
print("$P4V" in meta)
print(meta.get("$P4N"))
print(meta.get("$P4V"))

print("\n=== Trying with different string types ===")
print(meta.get(b"$P4N"))
print(meta.get(b"$P4V"))