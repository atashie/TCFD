# ISIMIP Pipeline - Troubleshooting Log

This document tracks errors encountered during development and their solutions.

---

## Format

Each issue should be documented as:

```
### Issue Title
**Date**: YYYY-MM-DD
**Phase**: [Setup/Config/Search/Download/Processing/Visualization]
**Error Message**:
<error details>

**Cause**:
<root cause analysis>

**Solution**:
<how it was resolved>

**Prevention**:
<how to avoid in future>
```

---

## Known Issues & Solutions

### Array Indexing on 3D NumPy Arrays
**Date**: 2026-01-13
**Phase**: Testing
**Error Message**:
```
IndexError: index 3 is out of bounds for axis 1 with size 1
```

**Cause**: Test fixtures using `np.broadcast_to()` created read-only views with shape `(time, 1, 1)` instead of full-sized arrays. Direct indexing like `data[50:52, 3, 3]` failed.

**Solution**: Use `.copy()` after broadcasting to create writable full-sized array:
```python
base = 273.15 + 20 * np.sin(np.arange(len(time))[:, None, None] / 365)
data = np.broadcast_to(base, (len(time), len(lat), len(lon))).copy()
```

**Prevention**: Always use `.copy()` when broadcasting arrays that will be modified.

---

### Calendar Conversion with NumPy int64
**Date**: 2026-01-13
**Phase**: Processing
**Error Message**:
```
DTypePromotionError: The DType <class 'numpy.dtypes.DateTime64DType'> could not be promoted by <class 'numpy.dtypes.Int64DType'>
```

**Cause**: Calendar conversion in `alignment.py` tried to use numpy int64 values directly with `timedelta`, causing dtype incompatibility.

**Solution**: Convert numpy int64 to Python int before arithmetic:
```python
day_int = int(day_num)  # Convert numpy.int64 to Python int
new_date = ref_date + timedelta(days=day_int)
```

**Prevention**: Always convert numpy scalar types to Python native types when using with stdlib functions.

---

### Validation Dictionary Key Mismatch
**Date**: 2026-01-13
**Phase**: Testing
**Error Message**:
```
AssertionError: assert 'partially_missing_cells' in gaps
```

**Cause**: Test expected `"partially_missing_cells"` but `detect_spatial_gaps()` returned `"n_partially_missing_cells"` (with `n_` prefix).

**Solution**: Updated test assertion to match actual key name.

**Prevention**: Use consistent naming conventions across all validation functions.

---

### Template Entry
**Date**: 2026-01-12
**Phase**: Setup
**Error Message**:
```
(No errors yet - project initialization)
```

**Cause**: N/A

**Solution**: N/A

**Prevention**: N/A

---

## Common Issues Reference

### NetCDF/xarray Issues

#### Issue: "ValueError: cannot reindex or align along dimension"
**Typical Cause**: Mismatched coordinate values between datasets
**Solution**: Use `xr.align()` with `join='inner'` or interpolate to common grid

#### Issue: "MemoryError when loading large NetCDF"
**Typical Cause**: Loading entire file into memory
**Solution**: Use `chunks` parameter in `xr.open_dataset()` for dask-backed lazy loading

#### Issue: "KeyError: 'time'"
**Typical Cause**: Different variable naming conventions across models
**Solution**: Check actual dimension names with `ds.dims`, may be 'Time' or 'date'

### ISIMIP API Issues

#### Issue: "Connection timeout"
**Typical Cause**: Large query result sets or server load
**Solution**: Add pagination, increase timeout, or narrow query filters

#### Issue: "No results found"
**Typical Cause**: Incorrect specifier values or typos
**Solution**: Verify against `/api/v1/identifiers/` endpoint for valid values

### Download Issues

#### Issue: "SSL: CERTIFICATE_VERIFY_FAILED"
**Typical Cause**: Corporate proxy or outdated certificates
**Solution**: Update certifi package or configure custom CA bundle

#### Issue: "Partial download / corrupted file"
**Typical Cause**: Network interruption
**Solution**: Implement checksum verification and resume capability

### Processing Issues

#### Issue: "Theil-Sen slope returns NaN"
**Typical Cause**: Insufficient non-NaN data points
**Solution**: Check for minimum data requirements before computing

#### Issue: "Spearman correlation warning about ties"
**Typical Cause**: Many identical values in data
**Solution**: This is informational; results are still valid

### Visualization Issues

#### Issue: "Plotly figure too large for browser"
**Typical Cause**: Too many data points in interactive plot
**Solution**: Downsample data or use WebGL renderer

---

## Environment-Specific Notes

### Windows

- Use raw strings or forward slashes for paths: `r"C:\path"` or `"C:/path"`
- NetCDF4 may require separate HDF5 installation
- Consider using WSL for better compatibility

### Conda vs pip

- Prefer conda for netcdf4, xarray, dask (better binary dependency handling)
- Can mix: `conda install netcdf4 xarray` then `pip install isimip-client`

---

## Debug Commands

```bash
# Check NetCDF file structure
ncdump -h <file.nc>

# Verify xarray can read file
python -c "import xarray as xr; print(xr.open_dataset('<file.nc>'))"

# Test ISIMIP API connection
curl "https://data.isimip.org/api/v1/datasets/?simulation_round=ISIMIP3b&limit=1"

# Check you.com API
curl --request POST \
  --url https://api.you.com/v1/agents/runs \
  --header 'Authorization: Bearer <api-key>' \
  --header 'Content-Type: application/json' \
  --data '{"agent": "<agent-id>", "stream": false}'
```

---

## Performance Optimization Notes

### Large File Processing

1. Use chunked reading: `xr.open_mfdataset(..., chunks={'time': 100})`
2. Process in batches rather than loading all models at once
3. Use `dask.distributed` for parallel processing on multi-core systems

### Memory Management

1. Close datasets after use: `ds.close()` or use context managers
2. Delete large variables when no longer needed: `del large_array`
3. Monitor memory with `tracemalloc` module

---

## Links to External Resources

- [xarray documentation](https://docs.xarray.dev/)
- [ISIMIP data portal](https://data.isimip.org)
- [ISIMIP protocol documentation](https://protocol.isimip.org/)
- [NetCDF conventions](https://cfconventions.org/)
- [Plotly Python docs](https://plotly.com/python/)
