import pandas as pd

# Especifica la ruta del archivo Excel
ruta_archivo_excel = "C:/Users/simon/Documents/GitHub/Probe-Designer/sondas/NC_000004.12_20240324_213512/NC_000004.12.xlsx"

# Carga todas las hojas del archivo Excel en un diccionario de DataFrames
try:
    xls = pd.ExcelFile(ruta_archivo_excel)
    hojas_excel = {}
    for nombre_hoja in xls.sheet_names:
        hojas_excel[nombre_hoja] = pd.read_excel(ruta_archivo_excel, sheet_name=nombre_hoja)
    
    # Muestra el contenido de cada hoja
    for nombre_hoja, df in hojas_excel.items():
        print(f"Contenido de la hoja '{nombre_hoja}':")
        print(df)
        print("\n")  # Agrega un espacio en blanco entre las hojas
except FileNotFoundError:
    print(f"No se pudo encontrar el archivo '{ruta_archivo_excel}'")
except Exception as e:
    print(f"Ocurrió un error al cargar el archivo Excel: {e}")
