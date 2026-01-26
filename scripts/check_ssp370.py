"""Check for SSP370 availability in marine fisheries sector."""
from isimip_client.client import ISIMIPClient

client = ISIMIPClient()

# Search all marine fishery data
print("Checking all available scenarios in marine-fishery_global...")

response = client.datasets(
    simulation_round='ISIMIP3b',
    sector='marine-fishery_global',
)

# Paginate through all
all_datasets = response.get('results', [])
while response.get('next'):
    import re
    page_match = re.search(r'page=(\d+)', response['next'])
    if page_match:
        response = client.datasets(
            simulation_round='ISIMIP3b',
            sector='marine-fishery_global',
            page=int(page_match.group(1))
        )
        all_datasets.extend(response.get('results', []))
    else:
        break

# Collect all unique scenarios
scenarios = {}
for ds in all_datasets:
    scenario = ds.get('specifiers', {}).get('climate_scenario', 'unknown')
    if scenario not in scenarios:
        scenarios[scenario] = 0
    scenarios[scenario] += 1

print(f"\nTotal datasets: {len(all_datasets)}")
print("\nScenarios available:")
for s, count in sorted(scenarios.items()):
    print(f"  {s}: {count} datasets")

# Specifically check for ssp370
ssp370 = [ds for ds in all_datasets if 'ssp370' in ds.get('specifiers', {}).get('climate_scenario', '')]
print(f"\nSSP370 datasets found: {len(ssp370)}")
