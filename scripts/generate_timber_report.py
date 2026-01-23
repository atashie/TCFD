"""Generate interactive HTML report for timber location extraction data.

Creates a single-page application with:
- Global map with points and polygons (color-coded by selectable metric)
- Hazard filter toggles to show/hide overlapping hazards
- Tabs to switch between Decade view and Scenario view
- Time series comparison (select up to 4 locations)
- Redesigned data table with scenario toggle and trend coloring
- Modern dark-mode aesthetic
- Lazy loading via external JSON data file
"""

import csv
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def log(msg: str):
    """Print with timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_wkt_polygon(wkt: str) -> List[List[float]]:
    """Parse WKT POLYGON string to list of [lon, lat] coordinates.

    Args:
        wkt: WKT string like "POLYGON ((-121.75 40.85, -121.74 40.86, ...))"

    Returns:
        List of [lon, lat] pairs
    """
    if not wkt or not wkt.strip():
        return []

    match = re.search(r'POLYGON\s*\(\((.+?)\)\)', wkt, re.IGNORECASE)
    if not match:
        return []

    coords_str = match.group(1)
    coords = []
    for pair in coords_str.split(','):
        parts = pair.strip().split()
        if len(parts) >= 2:
            lon = float(parts[0])
            lat = float(parts[1])
            coords.append([lon, lat])

    return coords


def normalize_longitude(lon: float) -> float:
    """Normalize longitude to -180 to 180 range.

    Handles multiple wraps (e.g., -452 becomes -92).

    Args:
        lon: Longitude value (may be outside normal range)

    Returns:
        Normalized longitude in -180 to 180 range
    """
    while lon < -180:
        lon += 360
    while lon > 180:
        lon -= 360
    return lon


def normalize_polygon_coords(coords: List[List[float]]) -> List[List[float]]:
    """Normalize all polygon coordinates to -180 to 180 range."""
    return [[normalize_longitude(lon), lat] for lon, lat in coords]


def load_location_polygons(csv_path: Path) -> Dict[str, List[List[float]]]:
    """Load polygon coordinates from source locations CSV.

    Now loads both poly rows AND region rows (which have polygon data after extraction).

    Args:
        csv_path: Path to ex-locations-timber.csv

    Returns:
        Dictionary mapping location names (with suffix) to polygon coordinates
    """
    polygons = {}

    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            location_name = row.get('Location Name', '').strip()
            polygon_wkt = row.get('Polygon', '').strip()
            region = row.get('Region', '').strip()

            if polygon_wkt:
                coords = parse_wkt_polygon(polygon_wkt)
                if coords:
                    coords = normalize_polygon_coords(coords)
                    # Determine key - handle both cases:
                    # 1. Location name already has suffix (e.g., "Shasta...(poly)")
                    # 2. Location name has no suffix (e.g., "Louisiana Loblolly Pine")
                    if location_name.endswith('(point)') or location_name.endswith('(poly)') or location_name.endswith('(region)'):
                        # Already has suffix, use directly
                        key = location_name
                    elif region:
                        # Has region column = region row
                        key = f"{location_name} (region)"
                    else:
                        # Polygon but no region = poly row
                        key = f"{location_name} (poly)"
                    polygons[key] = coords

    return polygons


def load_extraction_data(csv_path: Path) -> Tuple[List[Dict], Dict]:
    """Load extraction data from CSV export.

    Args:
        csv_path: Path to timber_locations_*.csv

    Returns:
        Tuple of (rows as list of dicts, metadata dict)
    """
    rows = []
    locations = set()
    scenarios = set()
    decades = set()
    hazards = set()

    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for field in ['Lat', 'Lon', 'Raw_Hazard_Value', 'Percentile_Score',
                         'Decadal_Trend_Strength', 'Decadal_Trend_Significance',
                         'Long_Term_Trend_Strength', 'Long_Term_Trend_Significance',
                         'Relative_Hazard_Score_Number', 'Raw_Hazard_Value_25th',
                         'Raw_Hazard_Value_75th']:
                if row.get(field) and row[field] != 'NA':
                    try:
                        row[field] = float(row[field])
                    except ValueError:
                        row[field] = None
                else:
                    row[field] = None

            # Convert decade to int
            if row.get('Decade'):
                row['Decade'] = int(row['Decade'])

            rows.append(row)
            locations.add(row['Location'])
            scenarios.add(row['Scenario'])
            decades.add(row['Decade'])
            if row.get('Hazard_Measure'):
                hazards.add(row['Hazard_Measure'])

    metadata = {
        'locations': sorted(list(locations)),
        'scenarios': sorted(list(scenarios)),
        'decades': sorted(list(decades)),
        'hazards': sorted(list(hazards)),
        'hazard': rows[0]['Hazard'] if rows else 'Timber',
        'hazard_measure': rows[0]['Hazard_Measure'] if rows else '',
        'generated': datetime.now().isoformat()
    }

    return rows, metadata


def build_location_info(rows: List[Dict], polygons: Dict[str, List[List[float]]]) -> List[Dict]:
    """Build location info list with coordinates and types.

    Args:
        rows: Extraction data rows
        polygons: Dictionary of polygon coordinates

    Returns:
        List of location info dicts
    """
    locations = {}

    for row in rows:
        loc_name = row['Location']
        if loc_name in locations:
            continue

        lat = row.get('Lat')
        lon = row.get('Lon')
        hazard_measure = row.get('Hazard_Measure', '')

        # Determine location type from name suffix
        if '(point)' in loc_name:
            loc_type = 'point'
        elif '(poly)' in loc_name:
            loc_type = 'polygon'
        elif '(region)' in loc_name:
            loc_type = 'region'
        else:
            loc_type = 'point' if lat and lon else 'unknown'

        loc_info = {
            'name': loc_name,
            'type': loc_type,
            'lat': lat,
            'lon': lon,
            'hazard': hazard_measure
        }

        # Add polygon coordinates if available
        if loc_name in polygons:
            loc_info['polygonCoords'] = polygons[loc_name]
            coords = polygons[loc_name]
            if coords:
                loc_info['centroidLon'] = sum(c[0] for c in coords) / len(coords)
                loc_info['centroidLat'] = sum(c[1] for c in coords) / len(coords)
                # Use centroid for region/polygon without point coords
                if loc_info['lat'] is None:
                    loc_info['lat'] = loc_info['centroidLat']
                if loc_info['lon'] is None:
                    loc_info['lon'] = loc_info['centroidLon']

        locations[loc_name] = loc_info

    return list(locations.values())


def generate_html(
    rows: List[Dict],
    metadata: Dict,
    locations: List[Dict],
    output_path: Path
):
    """Generate the HTML report file and data.json.

    Args:
        rows: Extraction data rows
        metadata: Metadata dictionary
        locations: Location info list
        output_path: Output HTML file path
    """
    # Write data to separate JSON file (for HTTP serving)
    data_json_path = output_path.parent / 'data.json'
    data_obj = {
        'rows': rows,
        'locations': locations,
        'metadata': metadata
    }
    with open(data_json_path, 'w', encoding='utf-8') as f:
        json.dump(data_obj, f)
    log(f"Wrote data to: {data_json_path}")

    # Also create embedded version of data for inline script
    data_json_inline = json.dumps(data_obj)


    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Timber Location Extraction Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --bg-primary: #1a1a2e;
            --bg-secondary: #16213e;
            --bg-card: #0f3460;
            --bg-hover: #1f4287;
            --text-primary: #eaeaea;
            --text-secondary: #94a3b8;
            --accent: #3498db;
            --accent-hover: #2980b9;
            --border: #2a4a7a;
            --success: #27ae60;
            --warning: #f39c12;
            --danger: #e74c3c;
            --scenario-low: #27ae60;
            --scenario-mid: #f39c12;
            --scenario-high: #e74c3c;
        }}

        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
            min-height: 100vh;
        }}

        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            padding: 24px 32px;
            border-bottom: 1px solid var(--border);
        }}

        .header h1 {{
            font-size: 1.8rem;
            font-weight: 600;
            margin-bottom: 4px;
        }}

        .header .subtitle {{
            color: var(--text-secondary);
            font-size: 0.95rem;
        }}

        .controls {{
            display: flex;
            gap: 24px;
            padding: 16px 32px;
            background: var(--bg-secondary);
            border-bottom: 1px solid var(--border);
            flex-wrap: wrap;
            align-items: center;
        }}

        .view-toggle {{
            display: flex;
            gap: 0;
            border-radius: 6px;
            overflow: hidden;
            border: 1px solid var(--border);
        }}

        .view-toggle button {{
            padding: 8px 16px;
            border: none;
            background: var(--bg-card);
            color: var(--text-secondary);
            cursor: pointer;
            font-size: 0.9rem;
            transition: all 0.2s;
        }}

        .view-toggle button:hover {{
            background: var(--bg-hover);
        }}

        .view-toggle button.active {{
            background: var(--accent);
            color: white;
        }}

        .metric-selector {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}

        .metric-selector label {{
            color: var(--text-secondary);
            font-size: 0.9rem;
        }}

        .metric-selector select {{
            padding: 8px 12px;
            border: 1px solid var(--border);
            border-radius: 6px;
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 0.9rem;
            cursor: pointer;
        }}

        .metric-selector select:focus {{
            outline: none;
            border-color: var(--accent);
        }}


        .tab-bar {{
            display: flex;
            gap: 0;
            padding: 0 32px;
            background: var(--bg-secondary);
            border-bottom: 1px solid var(--border);
            overflow-x: auto;
        }}

        .tab-bar button {{
            padding: 12px 20px;
            border: none;
            border-bottom: 3px solid transparent;
            background: transparent;
            color: var(--text-secondary);
            cursor: pointer;
            font-size: 0.9rem;
            transition: all 0.2s;
            white-space: nowrap;
        }}

        .tab-bar button:hover {{
            color: var(--text-primary);
            background: rgba(255,255,255,0.05);
        }}

        .tab-bar button.active {{
            color: var(--accent);
            border-bottom-color: var(--accent);
        }}

        .content {{
            padding: 24px 32px;
        }}

        .section {{
            margin-bottom: 32px;
        }}

        .section h2 {{
            font-size: 1.1rem;
            font-weight: 500;
            margin-bottom: 16px;
            color: var(--text-secondary);
        }}

        .card {{
            background: var(--bg-card);
            border-radius: 8px;
            border: 1px solid var(--border);
            overflow: hidden;
        }}

        #map-container {{
            height: 500px;
        }}

        #timeseries-container {{
            height: 350px;
        }}

        .selection-info {{
            padding: 12px 16px;
            background: var(--bg-secondary);
            font-size: 0.85rem;
            color: var(--text-secondary);
            display: flex;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
        }}

        .selection-badge {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 4px 10px;
            background: var(--bg-hover);
            border-radius: 16px;
            font-size: 0.8rem;
        }}

        .selection-badge .remove {{
            cursor: pointer;
            opacity: 0.7;
            font-weight: bold;
        }}

        .selection-badge .remove:hover {{
            opacity: 1;
        }}

        .selection-badge.color-0 {{ border-left: 3px solid #3498db; }}
        .selection-badge.color-1 {{ border-left: 3px solid #e74c3c; }}
        .selection-badge.color-2 {{ border-left: 3px solid #27ae60; }}
        .selection-badge.color-3 {{ border-left: 3px solid #9b59b6; }}

        .table-controls {{
            padding: 12px 16px;
            background: var(--bg-secondary);
            border-bottom: 1px solid var(--border);
            display: flex;
            gap: 16px;
            align-items: center;
            flex-wrap: wrap;
        }}

        .table-controls input {{
            padding: 8px 12px;
            border: 1px solid var(--border);
            border-radius: 6px;
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 0.9rem;
            width: 200px;
        }}

        .table-controls input:focus {{
            outline: none;
            border-color: var(--accent);
        }}

        .scenario-toggle {{
            display: flex;
            gap: 0;
            border-radius: 6px;
            overflow: hidden;
            border: 1px solid var(--border);
        }}

        .scenario-toggle button {{
            padding: 6px 12px;
            border: none;
            background: var(--bg-card);
            color: var(--text-secondary);
            cursor: pointer;
            font-size: 0.85rem;
            transition: all 0.2s;
        }}

        .scenario-toggle button:hover {{
            background: var(--bg-hover);
        }}

        .scenario-toggle button.active {{
            color: white;
        }}

        .scenario-toggle button[data-scenario="Low Emissions"].active {{
            background: var(--scenario-low);
        }}

        .scenario-toggle button[data-scenario="Middle of the Road"].active {{
            background: var(--scenario-mid);
        }}

        .scenario-toggle button[data-scenario="High Emissions"].active {{
            background: var(--scenario-high);
        }}

        .table-wrapper {{
            overflow-x: auto;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.85rem;
        }}

        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}

        th {{
            background: var(--bg-secondary);
            font-weight: 500;
            color: var(--text-secondary);
            cursor: pointer;
            user-select: none;
            position: sticky;
            top: 0;
            z-index: 1;
        }}

        th:hover {{
            background: var(--bg-hover);
        }}

        th .sort-icon {{
            margin-left: 6px;
            opacity: 0.5;
        }}

        th.sorted .sort-icon {{
            opacity: 1;
        }}

        tr:hover {{
            background: rgba(255,255,255,0.03);
        }}

        td input[type="checkbox"] {{
            width: 16px;
            height: 16px;
            cursor: pointer;
        }}

        td input[type="checkbox"]:disabled {{
            cursor: not-allowed;
            opacity: 0.5;
        }}

        .trend-worse {{
            background: rgba(231, 76, 60, 0.3);
        }}

        .trend-better {{
            background: rgba(39, 174, 96, 0.3);
        }}

        .risk-1 {{ color: #27ae60; }}
        .risk-2 {{ color: #2ecc71; }}
        .risk-3 {{ color: #f39c12; }}
        .risk-4 {{ color: #e67e22; }}
        .risk-5 {{ color: #e74c3c; }}

        .footer {{
            padding: 16px 32px;
            background: var(--bg-secondary);
            border-top: 1px solid var(--border);
            font-size: 0.8rem;
            color: var(--text-secondary);
            text-align: center;
        }}

        .placeholder {{
            display: flex;
            align-items: center;
            justify-content: center;
            height: 100%;
            color: var(--text-secondary);
            font-size: 0.95rem;
        }}

        .map-section {{
            display: flex;
            gap: 16px;
        }}

        .location-list-panel {{
            width: 280px;
            flex-shrink: 0;
            background: var(--bg-secondary);
            border-radius: 8px;
            border: 1px solid var(--border);
            overflow: hidden;
            display: flex;
            flex-direction: column;
        }}

        .location-list-panel h3 {{
            padding: 12px 16px;
            font-size: 0.9rem;
            font-weight: 500;
            color: var(--text-secondary);
            border-bottom: 1px solid var(--border);
            margin: 0;
        }}

        .location-list-panel .list-wrapper {{
            flex: 1;
            overflow-y: auto;
            max-height: 460px;
        }}

        .location-list-panel table {{
            width: 100%;
            font-size: 0.8rem;
        }}

        .location-list-panel td {{
            padding: 8px 12px;
        }}

        .location-list-panel tr:hover {{
            background: rgba(255,255,255,0.05);
        }}

        .location-list-panel input[type="checkbox"] {{
            width: 14px;
            height: 14px;
            cursor: pointer;
        }}

        .timeseries-controls {{
            padding: 12px 16px;
            background: var(--bg-secondary);
            border-bottom: 1px solid var(--border);
            display: flex;
            gap: 24px;
            align-items: center;
            flex-wrap: wrap;
        }}

        .timeseries-controls label {{
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 0.85rem;
            color: var(--text-primary);
            cursor: pointer;
        }}

        .timeseries-controls .scenario-toggles {{
            display: flex;
            gap: 12px;
            margin-left: auto;
        }}

        .timeseries-controls input[type="checkbox"] {{
            width: 14px;
            height: 14px;
            cursor: pointer;
        }}

        .loading {{
            display: flex;
            align-items: center;
            justify-content: center;
            height: 100vh;
            font-size: 1.2rem;
        }}

        @media (max-width: 768px) {{
            .controls {{
                padding: 12px 16px;
            }}
            .content {{
                padding: 16px;
            }}
            .tab-bar {{
                padding: 0 16px;
            }}
        }}
    </style>
</head>
<body>
    <div id="loading" class="loading">Loading data...</div>
    <div id="app" style="display: none;">
        <header class="header">
            <h1>Timber Location Extraction Report</h1>
            <div class="subtitle" id="header-subtitle"></div>
        </header>

        <nav class="controls">
            <div class="view-toggle">
                <button data-view="decade" class="active">View by Decade</button>
                <button data-view="scenario">View by Scenario</button>
            </div>
            <div class="metric-selector">
                <label>Color by:</label>
                <select id="metric-dropdown">
                    <option value="Raw_Hazard_Value">Raw Hazard Value</option>
                    <option value="Percentile_Score">Percentile Score (1-100)</option>
                    <option value="Relative_Hazard_Score_Number">Relative Hazard (1-5)</option>
                    <option value="Decadal_Trend_Strength">Trend Strength</option>
                    <option value="Decadal_Trend_Significance">Trend Significance (p-value)</option>
                </select>
            </div>
        </nav>

        <nav class="tab-bar" id="tab-container"></nav>

        <main class="content">
            <section class="section">
                <div class="map-section">
                    <div class="location-list-panel">
                        <h3>Locations</h3>
                        <div class="list-wrapper">
                            <table id="location-list-table">
                                <tbody></tbody>
                            </table>
                        </div>
                    </div>
                    <div class="card" style="flex: 1;">
                        <div id="map-container"></div>
                    </div>
                </div>
            </section>

            <section class="section">
                <h2>Compare Locations (select from list above - max 4)</h2>
                <div class="card">
                    <div class="selection-info" id="selection-info">
                        <span>No locations selected</span>
                    </div>
                    <div class="timeseries-controls">
                        <label>
                            <input type="checkbox" id="toggle-uncertainty">
                            Show Uncertainty (Q25-Q75)
                        </label>
                        <span class="scenario-toggles">
                            <label>
                                <input type="checkbox" data-ts-scenario="Low Emissions">
                                Low
                            </label>
                            <label>
                                <input type="checkbox" data-ts-scenario="Middle of the Road" checked>
                                Mid
                            </label>
                            <label>
                                <input type="checkbox" data-ts-scenario="High Emissions">
                                High
                            </label>
                        </span>
                    </div>
                    <div id="timeseries-container"></div>
                </div>
            </section>

            <section class="section">
                <h2>Data Table</h2>
                <div class="card">
                    <div class="table-controls">
                        <input type="text" id="location-filter" placeholder="Filter by location...">
                        <div class="scenario-toggle" id="scenario-toggle"></div>
                    </div>
                    <div class="table-wrapper" style="max-height: 400px; overflow-y: auto;">
                        <table id="data-table">
                            <thead>
                                <tr>
                                    <th style="width: 40px;"></th>
                                    <th data-col="Location">Location <span class="sort-icon">↕</span></th>
                                    <th data-col="Hazard_Measure">Hazard <span class="sort-icon">↕</span></th>
                                    <th data-col="2020_raw">2020 (kg/m²) <span class="sort-icon">↕</span></th>
                                    <th data-col="2020_pct">2020 (%) <span class="sort-icon">↕</span></th>
                                    <th data-col="2030_raw">2030 (kg/m²) <span class="sort-icon">↕</span></th>
                                    <th data-col="2030_pct">2030 (%) <span class="sort-icon">↕</span></th>
                                    <th data-col="2040_raw">2040 (kg/m²) <span class="sort-icon">↕</span></th>
                                    <th data-col="2040_pct">2040 (%) <span class="sort-icon">↕</span></th>
                                    <th data-col="2050_raw">2050 (kg/m²) <span class="sort-icon">↕</span></th>
                                    <th data-col="2050_pct">2050 (%) <span class="sort-icon">↕</span></th>
                                    <th data-col="Long_Term_Trend">Trend (kg/m²/decade) <span class="sort-icon">↕</span></th>
                                </tr>
                            </thead>
                            <tbody></tbody>
                        </table>
                    </div>
                </div>
            </section>
        </main>

        <footer class="footer">
            Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Data source: ISIMIP climate projections
        </footer>
    </div>

    <script>
    // Embedded data (works with file:// URLs)
    const EMBEDDED_DATA = {data_json_inline};
    let DATA = null;

    // Application state
    const AppState = {{
        viewMode: 'decade',
        selectedTab: null,
        selectedMetric: 'Raw_Hazard_Value',
        selectedLocations: [],
        tableScenario: 'Low Emissions',
        sortColumn: null,
        sortDirection: 'asc',
        filterText: '',
        showUncertainty: false,
        enabledScenarios: new Set(['Middle of the Road']),

        subscribers: [],

        subscribe(fn) {{
            this.subscribers.push(fn);
        }},

        notify() {{
            this.subscribers.forEach(fn => fn(this));
        }},

        update(partial) {{
            Object.assign(this, partial);
            this.notify();
        }}
    }};

    // Color scales
    const COLOR_SCALES = {{
        Raw_Hazard_Value: 'Viridis',
        Percentile_Score: 'RdYlBu',
        Relative_Hazard_Score_Number: [[0, '#27ae60'], [0.25, '#2ecc71'], [0.5, '#f39c12'], [0.75, '#e67e22'], [1, '#e74c3c']],
        Decadal_Trend_Strength: 'RdBu',
        Decadal_Trend_Significance: 'YlOrRd_r'
    }};

    const SELECTION_COLORS = ['#3498db', '#e74c3c', '#27ae60', '#9b59b6'];
    const SCENARIO_STYLES = {{
        'Low Emissions': {{ dash: 'dash', width: 2 }},
        'Middle of the Road': {{ dash: 'solid', width: 2 }},
        'High Emissions': {{ dash: 'dot', width: 2 }}
    }};

    // Helper functions
    function filterData() {{
        let filtered = DATA.rows;

        if (AppState.viewMode === 'decade') {{
            filtered = filtered.filter(r => r.Decade == AppState.selectedTab);
        }} else {{
            filtered = filtered.filter(r => r.Scenario === AppState.selectedTab);
        }}

        if (AppState.filterText) {{
            const search = AppState.filterText.toLowerCase();
            filtered = filtered.filter(r => r.Location.toLowerCase().includes(search));
        }}

        return filtered;
    }}

    function getLocationData(compositeKey) {{
        // compositeKey format: "Location Name|Hazard_Measure"
        const [locationName, hazardMeasure] = compositeKey.split('|');
        return DATA.rows.filter(r =>
            r.Location === locationName && r.Hazard_Measure === hazardMeasure
        );
    }}

    function getLocationInfo(locationName) {{
        return DATA.locations.find(l => l.name === locationName);
    }}

    function deriveHigherIsBetter() {{
        // Determine if higher raw values are DESIRABLE (good for the asset/ecosystem)
        // The percentile scoring is RISK-CENTRIC: lower percentile = less risk = better
        //
        // For vegetation carbon: higher raw = more carbon = GOOD = lower percentile (less risk)
        // So if higher raw correlates with LOWER percentile, higher IS better
        //
        // For fire risk: higher raw = more fire = BAD = higher percentile (more risk)
        // So if higher raw correlates with HIGHER percentile, higher is worse
        const sample = DATA.rows.filter(r => r.Raw_Hazard_Value !== null && r.Percentile_Score !== null);
        if (sample.length < 2) return true; // default to higher is better

        const sorted = [...sample].sort((a, b) => a.Raw_Hazard_Value - b.Raw_Hazard_Value);
        const first = sorted[0];  // Lowest raw value
        const last = sorted[sorted.length - 1];  // Highest raw value

        // If higher raw = lower percentile (less risk), then higher IS better
        // If higher raw = higher percentile (more risk), then higher is worse
        return first.Percentile_Score > last.Percentile_Score;
    }}

    function formatValue(value, metric) {{
        if (value === null || value === undefined) return 'N/A';
        if (metric === 'Percentile_Score') return value.toFixed(1);
        if (metric === 'Relative_Hazard_Score_Number') return Math.round(value);
        if (metric && metric.includes('Trend')) return value.toExponential(2);
        return value.toFixed(2);
    }}

    function getShortHazardName(hazard) {{
        if (!hazard) return '';
        return hazard
            .replace('Vegetation Carbon (kg/m2)', 'Yield')
            .replace('Temperate Needleleaf Evergreen ', 'Temperate Needleleaf ')
            .replace('Temperate Broadleaf Summer Deciduous ', 'Temperate Broadleaf ');
    }}

    function toggleLocationSelection(locationName) {{
        const idx = AppState.selectedLocations.indexOf(locationName);
        if (idx >= 0) {{
            AppState.selectedLocations.splice(idx, 1);
        }} else if (AppState.selectedLocations.length < 4) {{
            AppState.selectedLocations.push(locationName);
        }}
        AppState.notify();
    }}

    // Render functions
    function renderScenarioToggle() {{
        const container = document.getElementById('scenario-toggle');
        container.innerHTML = DATA.metadata.scenarios.map(s => {{
            const shortName = s === 'Middle of the Road' ? 'Mid' : (s === 'Low Emissions' ? 'Low' : 'High');
            const active = s === AppState.tableScenario ? 'active' : '';
            return `<button data-scenario="${{s}}" class="${{active}}">${{shortName}}</button>`;
        }}).join('');

        container.querySelectorAll('button').forEach(btn => {{
            btn.addEventListener('click', () => {{
                container.querySelectorAll('button').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                AppState.tableScenario = btn.dataset.scenario;
                AppState.notify();
            }});
        }});
    }}

    function renderTabs() {{
        const container = document.getElementById('tab-container');
        let tabs;

        if (AppState.viewMode === 'decade') {{
            tabs = DATA.metadata.decades.map(d => ({{ id: String(d), label: d + 's' }}));
        }} else {{
            tabs = DATA.metadata.scenarios.map(s => ({{ id: s, label: s }}));
        }}

        container.innerHTML = tabs.map(t =>
            `<button data-tab="${{t.id}}" class="${{t.id === AppState.selectedTab ? 'active' : ''}}">${{t.label}}</button>`
        ).join('');

        container.querySelectorAll('button').forEach(btn => {{
            btn.addEventListener('click', () => {{
                AppState.update({{ selectedTab: btn.dataset.tab }});
            }});
        }});
    }}

    function renderMap() {{
        const filteredData = filterData();
        const traces = [];

        // Group data by location
        const locationGroups = {{}};
        filteredData.forEach(row => {{
            if (!locationGroups[row.Location]) {{
                locationGroups[row.Location] = [];
            }}
            locationGroups[row.Location].push(row);
        }});

        // Determine color range
        const values = filteredData.map(r => r[AppState.selectedMetric]).filter(v => v !== null);
        let cmin = Math.min(...values);
        let cmax = Math.max(...values);

        // Fixed ranges for specific metrics
        if (AppState.selectedMetric === 'Percentile_Score') {{
            cmin = 0;
            cmax = 100;
        }} else if (AppState.selectedMetric === 'Relative_Hazard_Score_Number') {{
            cmin = 1;
            cmax = 5;
        }} else if (AppState.selectedMetric.includes('Trend_Strength')) {{
            // Symmetric range for trend
            const absMax = Math.max(Math.abs(cmin), Math.abs(cmax));
            cmin = -absMax;
            cmax = absMax;
        }}

        // Render each location
        DATA.locations.forEach(loc => {{
            const locData = locationGroups[loc.name];
            if (!locData || locData.length === 0) return;

            const value = locData[0][AppState.selectedMetric];
            // Check if any composite key with this location is selected
            const isSelected = AppState.selectedLocations.some(key => key.startsWith(loc.name + '|'));
            const selIdx = AppState.selectedLocations.findIndex(key => key.startsWith(loc.name + '|'));

            // Render polygon if available (for poly AND region types)
            if (loc.polygonCoords && loc.polygonCoords.length > 0) {{
                const lons = loc.polygonCoords.map(c => c[0]);
                const lats = loc.polygonCoords.map(c => c[1]);
                // Close the polygon
                lons.push(lons[0]);
                lats.push(lats[0]);

                const normalizedValue = (value - cmin) / (cmax - cmin || 1);

                traces.push({{
                    type: 'scattergeo',
                    lon: lons,
                    lat: lats,
                    mode: 'lines',
                    fill: 'toself',
                    fillcolor: getColorFromScale(normalizedValue, AppState.selectedMetric, 0.4),
                    line: {{
                        color: isSelected ? SELECTION_COLORS[selIdx] : getColorFromScale(normalizedValue, AppState.selectedMetric, 0.9),
                        width: isSelected ? 4 : 2
                    }},
                    text: loc.name,
                    customdata: Array(lons.length).fill(loc.name),
                    hovertemplate: '<b>' + loc.name + '</b><br>' + AppState.selectedMetric.replace(/_/g, ' ') + ': ' + formatValue(value, AppState.selectedMetric) + '<extra></extra>',
                    showlegend: false
                }});
            }} else if (loc.lat !== null && loc.lon !== null) {{
                // Point marker
                traces.push({{
                    type: 'scattergeo',
                    lon: [loc.lon],
                    lat: [loc.lat],
                    mode: 'markers',
                    marker: {{
                        color: [value],
                        colorscale: COLOR_SCALES[AppState.selectedMetric],
                        cmin: cmin,
                        cmax: cmax,
                        size: isSelected ? 18 : 12,
                        line: {{
                            color: isSelected ? SELECTION_COLORS[selIdx] : 'rgba(255,255,255,0.5)',
                            width: isSelected ? 3 : 1
                        }},
                        showscale: false
                    }},
                    text: [loc.name],
                    customdata: [loc.name],
                    hovertemplate: '<b>%{{text}}</b><br>' + AppState.selectedMetric.replace(/_/g, ' ') + ': %{{marker.color:.3f}}<extra></extra>',
                    showlegend: false
                }});
            }}
        }});

        // Add colorbar trace
        traces.push({{
            type: 'scattergeo',
            lon: [null],
            lat: [null],
            mode: 'markers',
            marker: {{
                color: [cmin, cmax],
                colorscale: COLOR_SCALES[AppState.selectedMetric],
                cmin: cmin,
                cmax: cmax,
                size: 0,
                showscale: true,
                colorbar: {{
                    title: {{
                        text: AppState.selectedMetric.replace(/_/g, ' '),
                        font: {{ color: '#eaeaea' }}
                    }},
                    tickfont: {{ color: '#94a3b8' }},
                    bgcolor: 'rgba(0,0,0,0)'
                }}
            }},
            showlegend: false,
            hoverinfo: 'skip'
        }});

        const layout = {{
            geo: {{
                showland: true,
                landcolor: '#1a1a2e',
                showocean: true,
                oceancolor: '#16213e',
                showcoastlines: true,
                coastlinecolor: '#2a4a7a',
                showframe: false,
                projection: {{ type: 'natural earth' }},
                bgcolor: 'rgba(0,0,0,0)'
            }},
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: {{ t: 10, b: 10, l: 10, r: 50 }},
            dragmode: 'pan'
        }};

        Plotly.react('map-container', traces, layout, {{ responsive: true }});
    }}

    function renderLocationList() {{
        const tbody = document.querySelector('#location-list-table tbody');
        const filteredData = filterData();

        // Group by location+hazard composite key
        const locationHazardValues = {{}};
        filteredData.forEach(row => {{
            const key = row.Location + '|' + row.Hazard_Measure;
            if (!locationHazardValues[key]) {{
                locationHazardValues[key] = {{
                    location: row.Location,
                    hazard: row.Hazard_Measure,
                    value: row[AppState.selectedMetric]
                }};
            }}
        }});

        // Sort by location name then hazard
        const sortedEntries = Object.entries(locationHazardValues).sort((a, b) => {{
            const locCmp = a[1].location.localeCompare(b[1].location);
            if (locCmp !== 0) return locCmp;
            return (a[1].hazard || '').localeCompare(b[1].hazard || '');
        }});

        tbody.innerHTML = sortedEntries.map(([key, data]) => {{
            const isSelected = AppState.selectedLocations.includes(key);
            const isDisabled = !isSelected && AppState.selectedLocations.length >= 4;
            // Extract type from suffix (point/poly/region)
            const typeMatch = data.location.match(/\\((point|poly|region)\\)$/);
            const locType = typeMatch ? typeMatch[1] : '';
            const shortName = data.location.split(' (')[0];
            const shortHazard = getShortHazardName(data.hazard);
            const displayName = shortName + (locType ? ` (${{locType}})` : '') + (shortHazard ? ` [${{shortHazard}}]` : '');

            return `<tr data-location="${{key.replace(/"/g, '&quot;')}}">
                <td>
                    <input type="checkbox" class="location-checkbox"
                        ${{isSelected ? 'checked' : ''}}
                        ${{isDisabled ? 'disabled' : ''}}
                    >
                </td>
                <td title="${{data.location}} - ${{data.hazard}}">${{displayName}}</td>
                <td>${{formatValue(data.value, AppState.selectedMetric)}}</td>
            </tr>`;
        }}).join('');
    }}

    function getColorFromScale(normalizedValue, metric, alpha) {{
        const v = Math.max(0, Math.min(1, normalizedValue));
        let r, g, b;

        if (metric === 'Relative_Hazard_Score_Number') {{
            r = Math.round(39 + (231 - 39) * v);
            g = Math.round(174 - (174 - 76) * v);
            b = Math.round(96 - (96 - 60) * v);
        }} else if (metric.includes('Trend_Strength')) {{
            if (v < 0.5) {{
                const t = v * 2;
                r = Math.round(33 + (255 - 33) * t);
                g = Math.round(102 + (255 - 102) * t);
                b = Math.round(172 + (255 - 172) * t);
            }} else {{
                const t = (v - 0.5) * 2;
                r = Math.round(255 - (255 - 178) * t);
                g = Math.round(255 - (255 - 24) * t);
                b = Math.round(255 - (255 - 43) * t);
            }}
        }} else if (metric === 'Percentile_Score') {{
            if (v < 0.5) {{
                const t = v * 2;
                r = Math.round(49 + (255 - 49) * t);
                g = Math.round(130 + (255 - 130) * t);
                b = Math.round(189 - (189 - 0) * t);
            }} else {{
                const t = (v - 0.5) * 2;
                r = 255;
                g = Math.round(255 - (255 - 0) * t);
                b = Math.round(0 + (65 - 0) * t);
            }}
        }} else {{
            r = Math.round(68 + (253 - 68) * v);
            g = Math.round(1 + (231 - 1) * v);
            b = Math.round(84 + (37 - 84) * v);
        }}

        return `rgba(${{r}},${{g}},${{b}},${{alpha}})`;
    }}

    function renderTimeSeries() {{
        const container = document.getElementById('timeseries-container');
        const infoContainer = document.getElementById('selection-info');

        // Purge existing plot to prevent memory leaks from repeated renders
        try {{
            Plotly.purge('timeseries-container');
        }} catch (e) {{
            // Ignore errors if container doesn't have a plot yet
        }}

        if (AppState.selectedLocations.length === 0) {{
            container.innerHTML = '<div class="placeholder">Select locations from the list above to compare trends</div>';
            infoContainer.innerHTML = '<span>No locations selected - check boxes in the location list</span>';
            return;
        }}

        // Render selection badges (composite key format: "Location|Hazard")
        infoContainer.innerHTML = AppState.selectedLocations.map((compositeKey, idx) => {{
            const [locName, hazard] = compositeKey.split('|');
            const shortName = locName.split(' (')[0];
            const shortHazard = getShortHazardName(hazard);
            const badgeText = shortName + (shortHazard ? ` [${{shortHazard}}]` : '');
            return `<span class="selection-badge color-${{idx}}">
                ${{badgeText}}
                <span class="remove" onclick="toggleLocationSelection('${{compositeKey.replace(/'/g, "\\'")}}')">×</span>
            </span>`;
        }}).join('');

        const traces = [];

        AppState.selectedLocations.forEach((compositeKey, locIdx) => {{
            const locData = getLocationData(compositeKey);
            const [locName, hazard] = compositeKey.split('|');
            const color = SELECTION_COLORS[locIdx];
            const shortHazard = getShortHazardName(hazard);

            // Group by scenario - only show enabled scenarios
            DATA.metadata.scenarios.forEach((scenario, sIdx) => {{
                // Skip scenarios that are not enabled
                if (!AppState.enabledScenarios.has(scenario)) return;

                const scenarioData = locData
                    .filter(r => r.Scenario === scenario)
                    .sort((a, b) => a.Decade - b.Decade);

                if (scenarioData.length === 0) return;

                const style = SCENARIO_STYLES[scenario] || {{ dash: 'solid', width: 2 }};
                const legendName = locName.split(' (')[0] + (shortHazard ? ` [${{shortHazard}}]` : '') + ' - ' + scenario;

                traces.push({{
                    type: 'scatter',
                    x: scenarioData.map(r => r.Decade),
                    y: scenarioData.map(r => r[AppState.selectedMetric]),
                    mode: 'lines+markers',
                    name: legendName,
                    legendgroup: compositeKey,
                    line: {{
                        color: color,
                        dash: style.dash,
                        width: style.width
                    }},
                    marker: {{
                        color: color,
                        size: 6
                    }},
                    showlegend: sIdx === 0 || !traces.some(t => t.legendgroup === compositeKey)
                }});

                // Add 25th-75th percentile uncertainty bands for Raw_Hazard_Value
                // Only show if uncertainty toggle is enabled
                if (AppState.showUncertainty && AppState.selectedMetric === 'Raw_Hazard_Value') {{
                    const has25th = scenarioData.some(r => r.Raw_Hazard_Value_25th !== null);
                    const has75th = scenarioData.some(r => r.Raw_Hazard_Value_75th !== null);

                    if (has25th && has75th) {{
                        // Convert hex color to rgba
                        const hexToRgba = (hex, alpha) => {{
                            const r = parseInt(hex.slice(1, 3), 16);
                            const g = parseInt(hex.slice(3, 5), 16);
                            const b = parseInt(hex.slice(5, 7), 16);
                            return `rgba(${{r}}, ${{g}}, ${{b}}, ${{alpha}})`;
                        }};

                        // Upper bound (75th percentile)
                        traces.push({{
                            type: 'scatter',
                            x: scenarioData.map(r => r.Decade),
                            y: scenarioData.map(r => r.Raw_Hazard_Value_75th),
                            mode: 'lines',
                            line: {{ color: 'transparent' }},
                            legendgroup: compositeKey + '-band',
                            showlegend: false,
                            hoverinfo: 'skip'
                        }});

                        // Lower bound (25th percentile) with fill to previous trace
                        traces.push({{
                            type: 'scatter',
                            x: scenarioData.map(r => r.Decade),
                            y: scenarioData.map(r => r.Raw_Hazard_Value_25th),
                            mode: 'lines',
                            fill: 'tonexty',
                            fillcolor: hexToRgba(color, 0.15),
                            line: {{ color: 'transparent' }},
                            legendgroup: compositeKey + '-band',
                            showlegend: false,
                            hoverinfo: 'skip'
                        }});
                    }}
                }}
            }});
        }});

        const layout = {{
            xaxis: {{
                title: {{ text: 'Decade', font: {{ color: '#94a3b8' }} }},
                tickfont: {{ color: '#94a3b8' }},
                gridcolor: '#2a4a7a',
                linecolor: '#2a4a7a'
            }},
            yaxis: {{
                title: {{ text: AppState.selectedMetric.replace(/_/g, ' '), font: {{ color: '#94a3b8' }} }},
                tickfont: {{ color: '#94a3b8' }},
                gridcolor: '#2a4a7a',
                linecolor: '#2a4a7a'
            }},
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: {{ t: 30, b: 60, l: 80, r: 30 }},
            legend: {{
                orientation: 'h',
                y: -0.2,
                font: {{ color: '#94a3b8' }}
            }},
            annotations: [{{
                text: 'Line styles: Dashed = Low Emissions, Solid = Middle of Road, Dotted = High Emissions',
                showarrow: false,
                x: 0.5,
                y: 1.05,
                xref: 'paper',
                yref: 'paper',
                font: {{ size: 11, color: '#94a3b8' }}
            }}]
        }};

        Plotly.react('timeseries-container', traces, layout, {{ responsive: true }});
    }}

    function renderTable() {{
        const tbody = document.querySelector('#data-table tbody');
        const scenario = AppState.tableScenario;

        // Filter by scenario
        let filtered = DATA.rows.filter(r => r.Scenario === scenario);

        // Apply location filter
        if (AppState.filterText) {{
            const search = AppState.filterText.toLowerCase();
            filtered = filtered.filter(r => r.Location.toLowerCase().includes(search));
        }}

        // Group by Location + Hazard_Measure (composite key)
        const groups = {{}};
        filtered.forEach(row => {{
            const key = row.Location + '|' + row.Hazard_Measure;
            if (!groups[key]) {{
                groups[key] = {{
                    compositeKey: key,
                    Location: row.Location,
                    Hazard_Measure: row.Hazard_Measure,
                    decades: {{}},
                    Long_Term_Trend_Strength: row.Long_Term_Trend_Strength,
                    Long_Term_Trend_Significance: row.Long_Term_Trend_Significance
                }};
            }}
            // Store both raw value and percentile for each decade
            groups[key].decades[row.Decade] = {{
                raw: row.Raw_Hazard_Value,
                percentile: row.Percentile_Score
            }};
            // Use the most recent long-term trend
            if (row.Decade === 2090 || row.Decade === Math.max(...DATA.metadata.decades)) {{
                groups[key].Long_Term_Trend_Strength = row.Long_Term_Trend_Strength;
                groups[key].Long_Term_Trend_Significance = row.Long_Term_Trend_Significance;
            }}
        }});

        let data = Object.values(groups);

        // Sort
        if (AppState.sortColumn) {{
            data = [...data].sort((a, b) => {{
                let aVal, bVal;
                // Handle decade columns (e.g., "2020_raw", "2030_pct")
                const decadeMatch = AppState.sortColumn.match(/^(\\d{{4}})_(raw|pct)$/);
                if (decadeMatch) {{
                    const decade = parseInt(decadeMatch[1]);
                    const field = decadeMatch[2] === 'raw' ? 'raw' : 'percentile';
                    aVal = a.decades[decade]?.[field];
                    bVal = b.decades[decade]?.[field];
                }} else if (AppState.sortColumn === 'Long_Term_Trend') {{
                    aVal = a.Long_Term_Trend_Strength;
                    bVal = b.Long_Term_Trend_Strength;
                }} else {{
                    aVal = a[AppState.sortColumn];
                    bVal = b[AppState.sortColumn];
                }}
                if (aVal === null || aVal === undefined) aVal = AppState.sortDirection === 'asc' ? Infinity : -Infinity;
                if (bVal === null || bVal === undefined) bVal = AppState.sortDirection === 'asc' ? Infinity : -Infinity;
                if (typeof aVal === 'string') {{
                    return AppState.sortDirection === 'asc'
                        ? aVal.localeCompare(bVal)
                        : bVal.localeCompare(aVal);
                }}
                return AppState.sortDirection === 'asc' ? aVal - bVal : bVal - aVal;
            }});
        }}

        tbody.innerHTML = data.map(row => {{
            // Use composite key for selection tracking
            const isSelected = AppState.selectedLocations.includes(row.compositeKey);
            const isDisabled = !isSelected && AppState.selectedLocations.length >= 4;

            // Determine trend class based on derived higher_is_better
            let trendClass = '';
            const trendStrength = row.Long_Term_Trend_Strength;
            const trendSig = row.Long_Term_Trend_Significance;
            const higherIsBetter = deriveHigherIsBetter();
            if (trendStrength !== null && trendSig !== null && trendSig < 0.05) {{
                // For higher_is_better metrics (like carbon), negative trend = worse
                // For lower_is_better metrics (like risk), negative trend = better
                if (trendStrength < 0) {{
                    trendClass = higherIsBetter ? 'trend-worse' : 'trend-better';
                }} else {{
                    trendClass = higherIsBetter ? 'trend-better' : 'trend-worse';
                }}
            }}

            const shortHazard = getShortHazardName(row.Hazard_Measure);

            return `<tr>
                <td>
                    <input type="checkbox"
                        ${{isSelected ? 'checked' : ''}}
                        ${{isDisabled ? 'disabled' : ''}}
                        onchange="toggleLocationSelection('${{row.compositeKey.replace(/'/g, "\\'")}}')"
                    >
                </td>
                <td>${{row.Location}}</td>
                <td title="${{row.Hazard_Measure}}">${{shortHazard}}</td>
                <td>${{formatValue(row.decades[2020]?.raw, 'Raw_Hazard_Value')}}</td>
                <td>${{formatValue(row.decades[2020]?.percentile, 'Percentile_Score')}}</td>
                <td>${{formatValue(row.decades[2030]?.raw, 'Raw_Hazard_Value')}}</td>
                <td>${{formatValue(row.decades[2030]?.percentile, 'Percentile_Score')}}</td>
                <td>${{formatValue(row.decades[2040]?.raw, 'Raw_Hazard_Value')}}</td>
                <td>${{formatValue(row.decades[2040]?.percentile, 'Percentile_Score')}}</td>
                <td>${{formatValue(row.decades[2050]?.raw, 'Raw_Hazard_Value')}}</td>
                <td>${{formatValue(row.decades[2050]?.percentile, 'Percentile_Score')}}</td>
                <td class="${{trendClass}}">${{formatValue(trendStrength, 'trend')}}</td>
            </tr>`;
        }}).join('');

        // Update header sort indicators
        document.querySelectorAll('#data-table th[data-col]').forEach(th => {{
            th.classList.toggle('sorted', th.dataset.col === AppState.sortColumn);
            const icon = th.querySelector('.sort-icon');
            if (th.dataset.col === AppState.sortColumn) {{
                icon.textContent = AppState.sortDirection === 'asc' ? '↑' : '↓';
            }} else {{
                icon.textContent = '↕';
            }}
        }});
    }}

    // Initialize application
    async function initApp() {{
        try {{
            // Use embedded data (works with file:// URLs)
            DATA = EMBEDDED_DATA;

            // Initialize state
            AppState.selectedTab = String(DATA.metadata.decades[0]);
            AppState.tableScenario = DATA.metadata.scenarios[0] || 'Low Emissions';

            // Update header
            document.getElementById('header-subtitle').textContent =
                `${{DATA.metadata.locations.length}} locations | ${{DATA.metadata.scenarios.length}} scenarios | ${{Math.min(...DATA.metadata.decades)}}s-${{Math.max(...DATA.metadata.decades)}}s`;

            // Hide loading, show app
            document.getElementById('loading').style.display = 'none';
            document.getElementById('app').style.display = 'block';

            // Setup event handlers
            document.querySelectorAll('.view-toggle button').forEach(btn => {{
                btn.addEventListener('click', () => {{
                    const newMode = btn.dataset.view;
                    const newTab = newMode === 'decade'
                        ? String(DATA.metadata.decades[0])
                        : DATA.metadata.scenarios[0];

                    document.querySelectorAll('.view-toggle button').forEach(b => b.classList.remove('active'));
                    btn.classList.add('active');

                    AppState.update({{ viewMode: newMode, selectedTab: newTab }});
                }});
            }});

            document.getElementById('metric-dropdown').addEventListener('change', (e) => {{
                AppState.update({{ selectedMetric: e.target.value }});
            }});

            document.getElementById('location-filter').addEventListener('input', (e) => {{
                AppState.filterText = e.target.value;
                AppState.notify();
            }});

            document.querySelectorAll('#data-table th[data-col]').forEach(th => {{
                th.addEventListener('click', () => {{
                    const col = th.dataset.col;
                    if (AppState.sortColumn === col) {{
                        AppState.sortDirection = AppState.sortDirection === 'asc' ? 'desc' : 'asc';
                    }} else {{
                        AppState.sortColumn = col;
                        AppState.sortDirection = 'asc';
                    }}
                    AppState.notify();
                }});
            }});

            // Time series toggles - uncertainty
            document.getElementById('toggle-uncertainty').addEventListener('change', (e) => {{
                AppState.showUncertainty = e.target.checked;
                AppState.notify();
            }});

            // Time series toggles - scenarios
            document.querySelectorAll('[data-ts-scenario]').forEach(cb => {{
                cb.addEventListener('change', () => {{
                    const scenario = cb.dataset.tsScenario;
                    if (cb.checked) {{
                        AppState.enabledScenarios.add(scenario);
                    }} else {{
                        AppState.enabledScenarios.delete(scenario);
                    }}
                    AppState.notify();
                }});
            }});

            // Location list checkbox event delegation (more robust than inline handlers)
            document.getElementById('location-list-table').addEventListener('change', (e) => {{
                if (e.target.classList.contains('location-checkbox')) {{
                    const loc = e.target.closest('tr').dataset.location;
                    if (loc) {{
                        toggleLocationSelection(loc);
                    }}
                }}
            }});

            // Subscribe to state changes
            AppState.subscribe(() => {{
                renderTabs();
                renderMap();
                renderLocationList();
                renderTimeSeries();
                renderTable();
            }});

            // Initial render
            renderScenarioToggle();
            renderTabs();
            renderMap();
            renderLocationList();
            renderTimeSeries();
            renderTable();

        }} catch (error) {{
            document.getElementById('loading').textContent = 'Error loading data: ' + error.message;
            console.error('Failed to load data:', error);
        }}
    }}

    // Start the app
    initApp();
    </script>
</body>
</html>'''

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    log(f"Generated report: {output_path}")


def main():
    """Main entry point."""
    # Paths
    project_root = Path(__file__).parent.parent
    extraction_csv = project_root / 'data' / 'exports' / 'timber_locations_20260123.csv'
    locations_csv = project_root / 'location-analyses' / 'ex-locations-timber.csv'
    output_dir = project_root / 'reports' / 'timber-extraction'
    output_path = output_dir / 'index.html'

    log("Loading data...")

    # Load polygon coordinates from source
    if locations_csv.exists():
        polygons = load_location_polygons(locations_csv)
        log(f"Loaded {len(polygons)} polygon definitions")
    else:
        log(f"Warning: Source locations file not found: {locations_csv}")
        polygons = {}

    # Load extraction data
    if not extraction_csv.exists():
        log(f"Error: Extraction CSV not found: {extraction_csv}")
        return

    rows, metadata = load_extraction_data(extraction_csv)
    log(f"Loaded {len(rows)} extraction rows")
    log(f"Locations: {metadata['locations']}")
    log(f"Scenarios: {metadata['scenarios']}")
    log(f"Decades: {metadata['decades']}")
    log(f"Hazards: {metadata['hazards']}")

    # Build location info
    locations = build_location_info(rows, polygons)
    log(f"Built info for {len(locations)} locations")

    # Generate HTML
    generate_html(rows, metadata, locations, output_path)

    log("Done!")


if __name__ == '__main__':
    main()
