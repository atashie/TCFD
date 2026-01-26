"""Search ISIMIP for large fish biomass data."""
from isimip_client.client import ISIMIPClient

client = ISIMIPClient()

# Search for marine fishery data
print("Searching ISIMIP3b marine-fishery_global...")
response = client.datasets(
    simulation_round='ISIMIP3b',
    sector='marine-fishery_global',
)

# Response is a dict with pagination
total_count = response.get('count', 0)
datasets = response.get('results', [])

print(f"\nFound {total_count} total datasets in marine-fishery_global sector")
print(f"Retrieved {len(datasets)} datasets in this page")

# Extract unique values
variables = set()
models = set()
scenarios = set()
gcms = set()

for ds in datasets:
    spec = ds.get('specifiers', {})
    if spec.get('variable'):
        variables.add(spec['variable'])
    if spec.get('model'):
        models.add(spec['model'])
    if spec.get('climate_scenario'):
        scenarios.add(spec['climate_scenario'])
    if spec.get('climate_forcing'):
        gcms.add(spec['climate_forcing'])

print(f"\n=== Available Data ===")
print(f"\nVariables: {sorted(variables)}")
print(f"\nModels: {sorted(models)}")
print(f"\nScenarios: {sorted(scenarios)}")
print(f"\nGCMs: {sorted(gcms)}")

# Look for size-class specific data
print("\n\n=== Size-Class Variables for Large Fish ===")
size_vars = [v for v in variables if 'log10' in v or 'size' in v.lower()]
print(f"Size-class variables: {size_vars}")

# Sample files for TCB
print("\n\n=== Sample TCB Datasets ===")
tcb_datasets = [ds for ds in datasets if ds.get('specifiers', {}).get('variable') == 'tcb']
print(f"Found {len(tcb_datasets)} TCB datasets")
for ds in tcb_datasets[:3]:
    spec = ds.get('specifiers', {})
    print(f"  - {spec.get('model')} / {spec.get('climate_forcing')} / {spec.get('climate_scenario')}")
    files = ds.get('files', [])
    if files:
        print(f"    File: {files[0].get('name')}")
        print(f"    Size: {files[0].get('size', 0) / 1e6:.1f} MB")

# Sample files for tclog10 (size-structured catch)
print("\n\n=== Sample tclog10 Datasets (Size-Structured Catch) ===")
tclog10_datasets = [ds for ds in datasets if ds.get('specifiers', {}).get('variable') == 'tclog10']
print(f"Found {len(tclog10_datasets)} tclog10 datasets")
for ds in tclog10_datasets[:3]:
    spec = ds.get('specifiers', {})
    print(f"  - {spec.get('model')} / {spec.get('climate_forcing')} / {spec.get('climate_scenario')}")
    files = ds.get('files', [])
    if files:
        print(f"    File: {files[0].get('name')}")
        print(f"    Size: {files[0].get('size', 0) / 1e6:.1f} MB")

# Get download URLs for recommended subset
print("\n\n=== Recommended Download: Large Fish Biomass ===")
print("Variable: tcb (Total Consumer Biomass) or tclog10 (size-structured)")
print("Scenarios: ssp126, ssp370, ssp585")
print("Models: boats, apecosm, dbpm, ecoocean, macroecological, zoomss")

# Filter for ssp scenarios (projection data)
ssp_tcb = [ds for ds in datasets
           if ds.get('specifiers', {}).get('variable') == 'tcb'
           and ds.get('specifiers', {}).get('climate_scenario', '').startswith('ssp')]
print(f"\nSSP projection TCB datasets: {len(ssp_tcb)}")

# Get unique model+gcm+scenario combinations for TCB
combos = set()
for ds in ssp_tcb:
    spec = ds.get('specifiers', {})
    combos.add((spec.get('model'), spec.get('climate_forcing'), spec.get('climate_scenario')))

print("\nAvailable combinations (model, GCM, scenario):")
for combo in sorted(combos):
    print(f"  {combo}")
