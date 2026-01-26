"""Inspect downloaded fish TCB data."""
import xarray as xr
from pathlib import Path

data_dir = Path("data/raw/fish-tcb")

# List all files
files = sorted(data_dir.glob("*.nc"))
print(f"Downloaded {len(files)} files:")
for f in files:
    print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")

# Inspect one file
print("\n\n=== Sample File Structure ===")
sample_file = files[0]
print(f"File: {sample_file.name}")

ds = xr.open_dataset(sample_file)
print(f"\nDimensions: {dict(ds.dims)}")
print(f"\nCoordinates:")
for coord in ds.coords:
    c = ds.coords[coord]
    print(f"  {coord}: {c.shape} ({c.dtype})")
    if c.size < 20:
        print(f"    values: {c.values}")
    else:
        print(f"    range: {c.values.min()} to {c.values.max()}")

print(f"\nData variables:")
for var in ds.data_vars:
    v = ds[var]
    print(f"  {var}: {v.dims} ({v.dtype})")
    if hasattr(v, 'long_name'):
        print(f"    long_name: {v.long_name}")
    if hasattr(v, 'units'):
        print(f"    units: {v.units}")

print(f"\nGlobal attributes:")
for attr in list(ds.attrs)[:10]:
    print(f"  {attr}: {ds.attrs[attr]}")

# Check TCB values
print("\n\n=== TCB Variable Summary ===")
tcb = ds['tcb']
print(f"Shape: {tcb.shape}")
print(f"Units: {tcb.attrs.get('units', 'N/A')}")
print(f"Long name: {tcb.attrs.get('long_name', 'N/A')}")

# Sample statistics (first time step, avoiding NaN)
print(f"\nValue range (first timestep):")
first_slice = tcb.isel(time=0).values
import numpy as np
valid = first_slice[~np.isnan(first_slice)]
print(f"  Min: {valid.min():.6f}")
print(f"  Max: {valid.max():.6f}")
print(f"  Mean: {valid.mean():.6f}")
print(f"  Valid cells: {len(valid)} / {first_slice.size}")

ds.close()
